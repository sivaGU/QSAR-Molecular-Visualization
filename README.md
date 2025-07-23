# QSAR Molecular Visualization Tool

A comprehensive web-based application for interactive 3D visualization and analysis of ERα and ERβ receptor-PFAS ligand structures with advanced QSAR/QSPR modeling capabilities.

## 🚀 Features

### Core Visualization
- **Interactive 3D Molecular Viewer**: Rotate, zoom, and explore molecular structures using NGL Viewer
- **ERα and ERβ Receptor Structures**: View 106 ligand complexes for each receptor type
- **Real-time Structure Loading**: Dynamic PDB file loading with protein-ligand visualization

### Data Analysis Dashboard
- **Statistical Summaries**: Comprehensive statistics for Alpha CE, Beta CE, Alpha T50, and Beta T50 datasets
- **Interactive Charts**: Bar charts, scatter plots, and distribution visualizations
- **Descriptor Analysis**: LogP, Molecular Weight, and PSA comparisons across datasets

### CE Ligand Comparison
- **Docking Score Analysis**: Compare docking scores between Alpha and Beta receptors
- **Descriptor Comparison**: Visualize chemical property differences
- **Performance Metrics**: Statistical analysis of commonly exposed ligands

### Chemical Descriptor Analysis
- **QSAR Coefficient Visualization**: Heatmaps showing descriptor importance
- **Top Influencer Charts**: Radar charts highlighting key molecular descriptors
- **Model Performance**: Interactive analysis of descriptor contributions

### QSAR Results
- **Model Performance Tracking**: R² values across refinement steps (Original, 10% Outliers Removed, 20% Outliers Removed)
- **Improvement Analysis**: Visual comparison of Alpha vs Beta model improvements
- **Statistical Insights**: Percentage improvements and key performance metrics

## 📁 Project Structure

```
Final GitHub Submission/
├── qsar_web_app.py              # Main Streamlit application
├── requirements.txt              # Python dependencies
├── run_app.bat                  # Windows batch file to run the app
├── run_app.ps1                  # PowerShell script to run the app
├── verify_setup.py              # Setup verification script
├── README.md                    # This file
├── INSTALLATION.md              # Detailed installation instructions
├── WHAT_YOU_GET.md              # Feature overview
├── .gitignore                   # Git ignore file
├── Alpha_CE_Combined/           # ERα commonly exposed ligand structures
├── Beta_CE_Combined/            # ERβ commonly exposed ligand structures
├── Alpha_T50_Combined/          # ERα T50 ligand structures
├── Beta_T50_Combined/           # ERβ T50 ligand structures
└── QSAR Submission/             # Supplementary data and tables
    ├── Supplementary Tables and Figures/
    └── Supplementary Codes/
```

## 🛠️ Installation

### Prerequisites
- Python 3.8 or higher
- pip package manager

### Quick Start
1. **Clone or download** this repository
2. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```
3. **Run the application**:
   ```bash
   streamlit run qsar_web_app.py
   ```
   
   Or use the provided scripts:
   - Windows: Double-click `run_app.bat`
   - PowerShell: Run `.\run_app.ps1`

### Alternative Installation Methods
See `INSTALLATION.md` for detailed installation instructions including:
- Virtual environment setup
- Docker installation
- Troubleshooting common issues

## 🎯 Usage

### Navigation
The application features a sidebar navigation with the following pages:
- **Home**: Overview and quick access to receptors
- **ERα Receptor**: Interactive 3D visualization of Alpha receptor structures
- **ERβ Receptor**: Interactive 3D visualization of Beta receptor structures
- **Data Analysis Dashboard**: Statistical analysis and visualizations
- **CE Ligand Comparison**: Comparison of commonly exposed ligands
- **Chemical Descriptor Analysis**: QSAR coefficient analysis
- **QSAR Results**: Model performance and improvement analysis
- **About**: Project information and documentation

### Key Features
1. **Molecular Visualization**: Select any of 106 ligands to view in 3D
2. **Data Analysis**: Explore statistical summaries and trends
3. **Model Performance**: Track QSAR model improvements
4. **Interactive Charts**: Hover, zoom, and explore data visualizations

## 📊 Data Sources

The application includes comprehensive datasets:
- **106 ERα ligand complexes** (Alpha_CE_Combined)
- **106 ERβ ligand complexes** (Beta_CE_Combined)
- **T50 datasets** for both receptors
- **QSAR model results** with refinement steps
- **Chemical descriptor coefficients** from manuscript analysis

## 🔬 Scientific Applications

This tool is designed for:
- **Drug Discovery**: Analyze ligand-receptor interactions
- **QSAR Modeling**: Evaluate model performance and improvements
- **Chemical Analysis**: Study molecular descriptors and properties
- **Educational Use**: Learn about molecular visualization and QSAR

## 🛡️ Technical Specifications

- **Framework**: Streamlit
- **3D Visualization**: NGL Viewer
- **Data Processing**: Pandas, NumPy
- **Visualization**: Plotly, Matplotlib
- **File Formats**: PDB, Excel, CSV

## 🤝 Contributing

For questions or contributions, please refer to the project documentation or contact the development team.


## 📄 License

This application is designed for research and educational use. Please cite the original analysis when using this data.

## 📞 Contact

For questions or issues: Contact Dr. Sivanesan Dakshanamurthy, PhD., MBA (sivanesan.dakshanamurthy@georgetown.edu)

