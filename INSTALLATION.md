# QSAR GUI Installation Guide

This guide provides step-by-step instructions for installing and running the QSAR Molecular Visualization Tool on different operating systems.

## Quick Start (All Platforms)

### Step 1: Download the Repository
1. Click the green "Code" button on GitHub
2. Select "Download ZIP"
3. Extract the ZIP file to a folder of your choice

### Step 2: Verify Installation
Run the verification script to check if everything is set up correctly:
```bash
python verify_setup.py
```

### Step 3: Launch the Application
Choose your preferred method based on your operating system.

---

## Windows Installation

### Method 1: One-Click Launch (Recommended)
1. **Download and extract** the repository
2. **Double-click** `run_app.bat`
3. **Wait** for dependencies to install
4. **Browser opens** automatically to `http://localhost:8501`

### Method 2: PowerShell
1. **Right-click** on `run_app.ps1`
2. **Select** "Run with PowerShell"
3. **Wait** for installation to complete
4. **Browser opens** automatically

### Method 3: Command Prompt
1. **Open Command Prompt** in the project folder
2. **Run**: `pip install -r requirements.txt`
3. **Run**: `streamlit run qsar_web_app.py`
4. **Open browser** to `http://localhost:8501`

### Troubleshooting Windows
- **Permission errors**: Run Command Prompt as Administrator
- **Python not found**: Install Python from [python.org](https://python.org)
- **Script execution policy**: Run `Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser` in PowerShell

---

## macOS Installation

### Prerequisites
1. **Install Python 3.9+**:
   ```bash
   # Using Homebrew (recommended)
   brew install python@3.11
   
   # Or download from python.org
   ```

2. **Install pip** (if not included):
   ```bash
   curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
   python3 get-pip.py
   ```

### Installation Steps
1. **Download and extract** the repository
2. **Open Terminal** in the project folder
3. **Install dependencies**:
   ```bash
   pip3 install -r requirements.txt
   ```
4. **Launch application**:
   ```bash
   streamlit run qsar_web_app.py
   ```
5. **Open browser** to `http://localhost:8501`

### Troubleshooting macOS
- **Permission errors**: Use `pip3` instead of `pip`
- **Python version**: Ensure you're using Python 3.9+
- **Browser issues**: Try Chrome or Safari

---

## Linux Installation

### Ubuntu/Debian
1. **Install Python 3.9+**:
   ```bash
   sudo apt update
   sudo apt install python3 python3-pip python3-venv
   ```

2. **Download and extract** the repository
3. **Open terminal** in the project folder
4. **Install dependencies**:
   ```bash
   pip3 install -r requirements.txt
   ```
5. **Launch application**:
   ```bash
   streamlit run qsar_web_app.py
   ```

### CentOS/RHEL/Fedora
1. **Install Python 3.9+**:
   ```bash
   # CentOS/RHEL
   sudo yum install python3 python3-pip
   
   # Fedora
   sudo dnf install python3 python3-pip
   ```

2. **Follow steps 2-5** from Ubuntu section above

### Troubleshooting Linux
- **Permission errors**: Use `pip3` and ensure proper permissions
- **Display issues**: Set `DISPLAY` environment variable if using remote access
- **Browser issues**: Install a modern browser (Chrome, Firefox)

---

## Advanced Installation

### Using Virtual Environment (Recommended)
```bash
# Create virtual environment
python -m venv qsar_env

# Activate virtual environment
# Windows:
qsar_env\Scripts\activate
# macOS/Linux:
source qsar_env/bin/activate

# Install dependencies
pip install -r requirements.txt

# Run application
streamlit run qsar_web_app.py
```

### Using Conda
```bash
# Create conda environment
conda create -n qsar_env python=3.11

# Activate environment
conda activate qsar_env

# Install dependencies
pip install -r requirements.txt

# Run application
streamlit run qsar_web_app.py
```

---

## Common Issues & Solutions

### Python Issues
- **"Python not found"**: Install Python 3.9+ from [python.org](https://python.org)
- **"pip not found"**: Install pip separately or use `python -m pip`
- **Version conflicts**: Use virtual environments

### Dependencies Issues
- **"Streamlit not found"**: Run `pip install -r requirements.txt`
- **Permission errors**: Use `pip install --user -r requirements.txt`
- **Network issues**: Check internet connection for package downloads

### Application Issues
- **"File not found"**: Ensure all data folders are present
- **"Port already in use"**: Close other applications or change port
- **Browser issues**: Try different browsers (Chrome recommended)

### Data Issues
- **Missing PDB files**: Re-download the repository
- **Corrupted files**: Check file integrity
- **Permission errors**: Ensure read permissions on data folders

---

## Getting Help

If you encounter issues:

1. **Run the verification script**: `python verify_setup.py`
2. **Check the troubleshooting section** in the main README
3. **Verify your Python version**: `python --version`
4. **Check your browser console** for JavaScript errors (F12)
5. **Ensure all 212 PDB files** are present in the data folders

For additional support, please include:
- Your operating system and version
- Python version (`python --version`)
- Error messages from the verification script
- Browser type and version

---

## Success!

Once you see the verification script pass all checks, you're ready to use the QSAR GUI! The application will provide interactive 3D visualization of 212 PFAS ligand-receptor complexes with full functionality. 
