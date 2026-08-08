#!/usr/bin/env python3
import sqlite3
import csv
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
CSV_PATH = BASE_DIR / "cell-count.csv"
DB_PATH = BASE_DIR / "cell-count.db"

CELL_COLUMNS = {
    "b_cell": "B cell",
    "cd8_t_cell": "CD8 T cell",
    "cd4_t_cell": "CD4 T cell",
    "nk_cell": "NK cell",
    "monocyte": "Monocyte",
}

EXPECTED_COLUMNS = {
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
    *CELL_COLUMNS.keys(),
}

SCHEMA = """
PRAGMA foreign_keys = ON;

CREATE TABLE projects (
    project_id INTEGER PRIMARY KEY,
    project_name TEXT NOT NULL UNIQUE
);

CREATE TABLE subjects (
    subject_id INTEGER PRIMARY KEY,
    subject_name TEXT NOT NULL UNIQUE,
    project_id INTEGER NOT NULL,
    condition TEXT NOT NULL,
    age INTEGER NOT NULL CHECK (age >= 0),
    sex TEXT NOT NULL,
    treatment TEXT NOT NULL,
    response TEXT,
    sample_type TEXT NOT NULL,
    FOREIGN KEY (project_id) REFERENCES projects(project_id)
);

CREATE TABLE samples (
    sample_id INTEGER PRIMARY KEY,
    sample_name TEXT NOT NULL UNIQUE,
    subject_id INTEGER NOT NULL,
    time_from_treatment_start INTEGER NOT NULL,
    FOREIGN KEY (subject_id) REFERENCES subjects(subject_id),
    UNIQUE (subject_id, time_from_treatment_start)
);

CREATE TABLE cell_types (
    cell_type_id INTEGER PRIMARY KEY,
    cell_type_key TEXT NOT NULL UNIQUE,
    cell_type_name TEXT NOT NULL UNIQUE
);

CREATE TABLE cell_counts (
    sample_id INTEGER NOT NULL,
    cell_type_id INTEGER NOT NULL,
    cell_count INTEGER NOT NULL CHECK (cell_count >= 0),
    PRIMARY KEY (sample_id, cell_type_id),
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id) ON DELETE CASCADE,
    FOREIGN KEY (cell_type_id) REFERENCES cell_types(cell_type_id)
);

CREATE INDEX idx_subjects_project_id ON subjects(project_id);
CREATE INDEX idx_samples_subject_id ON samples(subject_id);
CREATE INDEX idx_cell_counts_cell_type_id ON cell_counts(cell_type_id);
"""

def parse_int(value: str, column: str, row_number: int) -> int:
    """Parse a required integer and raise a useful error if it is invalid."""
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Row {row_number}: expected integer in {column!r}, got {value!r}"
        ) from exc


def nullable_text(value: str | None) -> str | None:
    """Convert a blank CSV value to SQL NULL."""
    if value is None:
        return None
    value = value.strip()
    return value if value else None


def load_database(csv_path: Path = CSV_PATH, db_path: Path = DB_PATH) -> None:
    if not csv_path.is_file():
        raise FileNotFoundError(
            f"Could not find {csv_path.name!r} next to {Path(__file__).name!r}."
        )

    # Recreate the database so every run produces a clean, deterministic result.
    if db_path.exists():
        db_path.unlink()

    connection = sqlite3.connect(db_path)
    connection.execute("PRAGMA foreign_keys = ON")

    try:
        with connection:
            connection.executescript(SCHEMA)

            connection.executemany(
                """
                INSERT INTO cell_types (cell_type_key, cell_type_name)
                VALUES (?, ?)
                """,
                CELL_COLUMNS.items(),
            )

            cell_type_ids = {
                key: cell_type_id
                for cell_type_id, key in connection.execute(
                    "SELECT cell_type_id, cell_type_key FROM cell_types"
                )
            }

            project_ids: dict[str, int] = {}
            subject_ids: dict[str, int] = {}

            with csv_path.open("r", encoding="utf-8-sig", newline="") as csv_file:
                reader = csv.DictReader(csv_file)

                if reader.fieldnames is None:
                    raise ValueError("CSV file has no header row.")

                actual_columns = set(reader.fieldnames)
                missing = EXPECTED_COLUMNS - actual_columns
                extra = actual_columns - EXPECTED_COLUMNS
                if missing or extra:
                    details = []
                    if missing:
                        details.append(f"missing columns: {sorted(missing)}")
                    if extra:
                        details.append(f"unexpected columns: {sorted(extra)}")
                    raise ValueError("CSV schema mismatch: " + "; ".join(details))

                for row_number, row in enumerate(reader, start=2):
                    project_name = row["project"].strip()
                    subject_name = row["subject"].strip()
                    sample_name = row["sample"].strip()

                    if project_name not in project_ids:
                        cursor = connection.execute(
                            "INSERT INTO projects (project_name) VALUES (?)",
                            (project_name,),
                        )
                        project_ids[project_name] = cursor.lastrowid

                    subject_metadata = (
                        project_ids[project_name],
                        row["condition"].strip(),
                        parse_int(row["age"], "age", row_number),
                        row["sex"].strip(),
                        row["treatment"].strip(),
                        nullable_text(row["response"]),
                        row["sample_type"].strip(),
                    )

                    if subject_name not in subject_ids:
                        cursor = connection.execute(
                            """
                            INSERT INTO subjects (
                                subject_name,
                                project_id,
                                condition,
                                age,
                                sex,
                                treatment,
                                response,
                                sample_type
                            )
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                            """,
                            (subject_name, *subject_metadata),
                        )
                        subject_ids[subject_name] = cursor.lastrowid
                    else:
                        # Protect against silently loading inconsistent
                        # subject-level metadata from later CSV rows.
                        stored = connection.execute(
                            """
                            SELECT
                                project_id,
                                condition,
                                age,
                                sex,
                                treatment,
                                response,
                                sample_type
                            FROM subjects
                            WHERE subject_id = ?
                            """,
                            (subject_ids[subject_name],),
                        ).fetchone()

                        if stored != subject_metadata:
                            raise ValueError(
                                f"Row {row_number}: inconsistent metadata "
                                f"for subject {subject_name!r}."
                            )

                    cursor = connection.execute(
                        """
                        INSERT INTO samples (
                            sample_name,
                            subject_id,
                            time_from_treatment_start
                        )
                        VALUES (?, ?, ?)
                        """,
                        (
                            sample_name,
                            subject_ids[subject_name],
                            parse_int(
                                row["time_from_treatment_start"],
                                "time_from_treatment_start",
                                row_number,
                            ),
                        ),
                    )
                    sample_id = cursor.lastrowid

                    counts = []
                    for column_name in CELL_COLUMNS:
                        count = parse_int(row[column_name], column_name, row_number)
                        if count < 0:
                            raise ValueError(
                                f"Row {row_number}: {column_name!r} cannot be negative."
                            )
                        counts.append(
                            (sample_id, cell_type_ids[column_name], count)
                        )

                    connection.executemany(
                        """
                        INSERT INTO cell_counts (
                            sample_id,
                            cell_type_id,
                            cell_count
                        )
                        VALUES (?, ?, ?)
                        """,
                        counts,
                    )

        sample_count = connection.execute(
            "SELECT COUNT(*) FROM samples"
        ).fetchone()[0]
        subject_count = connection.execute(
            "SELECT COUNT(*) FROM subjects"
        ).fetchone()[0]
        count_rows = connection.execute(
            "SELECT COUNT(*) FROM cell_counts"
        ).fetchone()[0]

        print(f"Created {db_path}")
        print(
            f"Loaded {sample_count:,} samples from {subject_count:,} subjects "
            f"and {count_rows:,} cell-count measurements."
        )

    except Exception:
        connection.close()
        if db_path.exists():
            db_path.unlink()
        raise
    else:
        connection.close()


if __name__ == "__main__":
    load_database()