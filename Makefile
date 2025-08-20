.PHONY: install run build build-linux build-linux-onefile build-win clean

VENV=.venv
PYTHON=$(VENV)/bin/python
PIP=$(VENV)/bin/pip

install:
	python3 -m venv $(VENV)
	. $(VENV)/bin/activate && $(PIP) install --upgrade pip
	. $(VENV)/bin/activate && $(PIP) install -r requirements.txt

run:
	. $(VENV)/bin/activate && $(PYTHON) magrf.py

# Legacy build (onefile, GUI). May fail on some distros due to libpython issues.
build:
	. $(VENV)/bin/activate && $(PIP) install pyinstaller pyinstaller-hooks-contrib
	. $(VENV)/bin/activate && pyinstaller --noconfirm --onefile --windowed \
		--collect-all matplotlib \
		--collect-all ttkbootstrap \
		 --hidden-import PIL._tkinter_finder \
		 --add-data "resources/icon.png:resources" \
		 --add-data "resources/icon.ico:resources" \
		 --icon resources/icon.png \
		--name calculadora-casamento magrf.py

# Recommended on Linux: onedir build (more robust w.r.t. shared libs)
build-linux:
	. $(VENV)/bin/activate && $(PIP) install pyinstaller pyinstaller-hooks-contrib
	. $(VENV)/bin/activate && pyinstaller --noconfirm --onedir --windowed \
		--collect-all matplotlib \
		--collect-all ttkbootstrap \
		 --hidden-import PIL._tkinter_finder \
		 --add-data "resources/icon.png:resources" \
		 --add-data "resources/icon.ico:resources" \
		 --icon resources/icon.png \
		--name calculadora-casamento magrf.py

# Try onefile on Linux again with extra collections
build-linux-onefile:
	. $(VENV)/bin/activate && $(PIP) install pyinstaller pyinstaller-hooks-contrib
	. $(VENV)/bin/activate && pyinstaller --noconfirm --onefile --windowed \
		--collect-all matplotlib \
		--collect-all ttkbootstrap \
		 --hidden-import PIL._tkinter_finder \
		 --add-data "resources/icon.png:resources" \
		 --add-data "resources/icon.ico:resources" \
		 --icon resources/icon.png \
		--name calculadora-casamento magrf.py

# Windows build must be done on Windows (or via CI). This target is a hint.
build-win:
	@echo "Build para Windows deve ser feito no Windows (ou via GitHub Actions)."
	@echo "Use o workflow em .github/workflows/build.yml ou rode:"
	@echo "  py -m venv .venv && .venv\\Scripts\\pip install -r requirements.txt pyinstaller"
	@echo "  .venv\\Scripts\\pyinstaller --noconfirm --onefile --windowed --collect-all matplotlib --collect-all ttkbootstrap --hidden-import PIL._tkinter_finder --add-data \"resources/icon.png;resources\" --add-data \"resources/icon.ico;resources\" --icon resources/icon.ico --name calculadora-casamento magrf.py"

clean:
	rm -rf build dist *.spec __pycache__
