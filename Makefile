.PHONY: setup pipeline dashboard

PYTHON ?= python3
PIPELINE_SCRIPT ?= load_data.py
DASHBOARD_APP ?= dashboard.py

setup:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -r requirements.txt

pipeline:
	$(PYTHON) $(PIPELINE_SCRIPT)

dashboard:
	$(PYTHON) -m streamlit run $(DASHBOARD_APP)
