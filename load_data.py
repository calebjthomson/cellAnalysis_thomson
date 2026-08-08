#!/usr/bin/env python3
import sqlite3
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu

BASE_DIR = Path(__file__).resolve().parent
CSV_PATH = BASE_DIR / "cell-count.csv"
DB_PATH = BASE_DIR / "cell-count.db"
BOXPLOT_PATH = BASE_DIR / "miraclib_response_boxplots.png"
STATS_PATH = BASE_DIR / "miraclib_response_statistics.csv"
BASELINE_TIME = 0

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

def get_frequency_summary(connection: sqlite3.Connection) -> pd.DataFrame:
    """Return relative cell-population frequencies for every sample."""
    query = """
        WITH sample_totals AS (
            SELECT
                sample_id,
                SUM(cell_count) AS total_count
            FROM cell_counts
            GROUP BY sample_id
        )
        SELECT
            s.sample_name AS sample,
            st.total_count,
            ct.cell_type_key AS population,
            cc.cell_count AS count,
            CASE
                WHEN st.total_count = 0 THEN 0.0
                ELSE 100.0 * cc.cell_count / st.total_count
            END AS percentage
        FROM cell_counts AS cc
        JOIN samples AS s
            ON s.sample_id = cc.sample_id
        JOIN cell_types AS ct
            ON ct.cell_type_id = cc.cell_type_id
        JOIN sample_totals AS st
            ON st.sample_id = cc.sample_id
        ORDER BY s.sample_name, ct.cell_type_id
    """
    return pd.read_sql_query(query, connection)


def display_frequency_summary(summary: pd.DataFrame) -> None:
    """Print the relative-frequency summary table."""
    print("\nRelative frequency summary")
    print(
        summary.to_string(
            index=False,
            formatters={"percentage": lambda value: f"{value:.2f}"},
        )
    )


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    """Return Benjamini-Hochberg FDR-adjusted p-values."""
    count = len(p_values)
    order = sorted(range(count), key=lambda index: p_values[index])
    adjusted = [0.0] * count
    running_min = 1.0

    for rank_index in range(count - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        candidate = p_values[original_index] * count / rank
        running_min = min(running_min, candidate)
        adjusted[original_index] = min(running_min, 1.0)

    return adjusted


def analyze_miraclib_response(
    connection: sqlite3.Connection,
    summary: pd.DataFrame,
) -> pd.DataFrame:
    """Compare baseline PBMC frequencies for miraclib responders/non-responders.

    Only baseline samples are used because the stated goal is prediction of
    treatment response. Including post-treatment measurements would leak
    information that is unavailable at prediction time.
    """
    metadata_query = """
        SELECT
            s.sample_name AS sample,
            s.time_from_treatment_start,
            sub.subject_name AS subject,
            sub.condition,
            sub.treatment,
            sub.response,
            sub.sample_type
        FROM samples AS s
        JOIN subjects AS sub
            ON sub.subject_id = s.subject_id
    """
    metadata = pd.read_sql_query(metadata_query, connection)
    analysis = summary.merge(metadata, on="sample", how="left")

    analysis = analysis[
        (analysis["condition"].str.lower() == "melanoma")
        & (analysis["treatment"].str.lower() == "miraclib")
        & (analysis["sample_type"].str.upper() == "PBMC")
        & (analysis["response"].str.lower().isin(["yes", "no"]))
        & (analysis["time_from_treatment_start"] == BASELINE_TIME)
    ].copy()

    if analysis.empty:
        raise ValueError(
            "No baseline PBMC samples were found for melanoma patients "
            "receiving miraclib."
        )

    response_counts = (
        analysis[["subject", "response"]]
        .drop_duplicates()
        ["response"]
        .value_counts()
    )
    print(
        "\nMiraclib response analysis: baseline melanoma PBMC samples "
        f"({BASELINE_TIME} days from treatment start)"
    )
    print(
        f"Responders: {response_counts.get('yes', 0)} subjects; "
        f"non-responders: {response_counts.get('no', 0)} subjects."
    )

    results = []
    raw_p_values = []

    for population in CELL_COLUMNS:
        population_data = analysis[analysis["population"] == population]
        responders = population_data.loc[
            population_data["response"].str.lower() == "yes", "percentage"
        ].astype(float)
        nonresponders = population_data.loc[
            population_data["response"].str.lower() == "no", "percentage"
        ].astype(float)

        if responders.empty or nonresponders.empty:
            raise ValueError(
                f"Both response groups are required for population {population!r}."
            )

        test = mannwhitneyu(
            responders,
            nonresponders,
            alternative="two-sided",
            method="auto",
        )

        # Rank-biserial correlation derived from Mann-Whitney U.
        # Positive values indicate higher frequencies among responders.
        effect_size = (
            2.0 * float(test.statistic) / (len(responders) * len(nonresponders))
            - 1.0
        )

        raw_p_values.append(float(test.pvalue))
        results.append(
            {
                "population": population,
                "n_responders": len(responders),
                "n_nonresponders": len(nonresponders),
                "responder_median_pct": responders.median(),
                "nonresponder_median_pct": nonresponders.median(),
                "median_difference_pct_points": (
                    responders.median() - nonresponders.median()
                ),
                "mann_whitney_u": float(test.statistic),
                "p_value": float(test.pvalue),
                "rank_biserial_correlation": effect_size,
            }
        )

    adjusted_p_values = benjamini_hochberg(raw_p_values)
    for row, adjusted_p in zip(results, adjusted_p_values):
        row["adjusted_p_value"] = adjusted_p
        row["significant_fdr_0_05"] = adjusted_p < 0.05

    stats = pd.DataFrame(results)
    stats.to_csv(STATS_PATH, index=False)

    print("\nResponder vs non-responder statistics")
    printable = stats.copy()
    numeric_columns = [
        "responder_median_pct",
        "nonresponder_median_pct",
        "median_difference_pct_points",
        "mann_whitney_u",
        "p_value",
        "adjusted_p_value",
        "rank_biserial_correlation",
    ]
    for column in numeric_columns:
        printable[column] = printable[column].map(lambda value: f"{value:.6g}")
    print(printable.to_string(index=False))

    significant = stats.loc[
        stats["significant_fdr_0_05"], "population"
    ].tolist()
    if significant:
        print(
            "\nSignificant populations after Benjamini-Hochberg correction "
            "(FDR < 0.05): "
            + ", ".join(significant)
        )
    else:
        print(
            "\nNo cell populations were significant after "
            "Benjamini-Hochberg correction (FDR < 0.05)."
        )

    plot_response_boxplots(analysis)
    print(f"Saved statistics to {STATS_PATH}")
    print(f"Saved boxplots to {BOXPLOT_PATH}")

    return stats


def plot_response_boxplots(analysis: pd.DataFrame) -> None:
    """Save responder/non-responder boxplots for all immune populations."""
    populations = list(CELL_COLUMNS)
    figure, axes = plt.subplots(
        1,
        len(populations),
        figsize=(16, 5),
        sharey=True,
    )

    for axis, population in zip(axes, populations):
        population_data = analysis[analysis["population"] == population]
        responders = population_data.loc[
            population_data["response"].str.lower() == "yes", "percentage"
        ].astype(float)
        nonresponders = population_data.loc[
            population_data["response"].str.lower() == "no", "percentage"
        ].astype(float)

        axis.boxplot(
            [responders, nonresponders],
            tick_labels=["Responder", "Non-responder"],
        )
        axis.set_title(CELL_COLUMNS[population])
        axis.tick_params(axis="x", labelrotation=30)

    axes[0].set_ylabel("Relative frequency (%)")
    figure.suptitle(
        "Baseline PBMC cell frequencies: miraclib responders vs non-responders"
    )
    figure.tight_layout()
    figure.savefig(BOXPLOT_PATH, dpi=200, bbox_inches="tight")
    plt.close(figure)



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

        summary = get_frequency_summary(connection)
        display_frequency_summary(summary)
        analyze_miraclib_response(connection, summary)

    except Exception:
        connection.close()
        if db_path.exists():
            db_path.unlink()
        raise
    else:
        connection.close()


if __name__ == "__main__":
    load_database()