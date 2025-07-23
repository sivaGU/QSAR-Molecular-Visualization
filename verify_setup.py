#!/usr/bin/env python3
"""
QSAR GUI Setup Verification Script
==================================

This script verifies that your QSAR GUI installation is working correctly.
Run this script to check if everything is set up properly before using the application.
"""

import sys
import subprocess
from pathlib import Path

def check_python_version():
    """Check if Python version is compatible."""
    print("🐍 Checking Python version...")
    version = sys.version_info
    if version.major == 3 and version.minor >= 9:
        print(f"✅ Python {version.major}.{version.minor}.{version.micro} - Compatible")
        return True
    else:
        print(f"❌ Python {version.major}.{version.minor}.{version.micro} - Requires Python 3.9+")
        return False

def check_dependencies():
    """Check if required dependencies are installed."""
    print("\n📦 Checking dependencies...")
    
    try:
        import streamlit
        print(f"✅ Streamlit {streamlit.__version__} - Installed")
        return True
    except ImportError:
        print("❌ Streamlit - Not installed")
        print("   Run: pip install -r requirements.txt")
        return False

def check_data_files():
    """Check if all data files are present."""
    print("\n📁 Checking data files...")
    
    expected_folders = ['Alpha_CE_Combined', 'Beta_CE_Combined', 'Alpha_T50_Combined', 'Beta_T50_Combined']
    total_files = 0
    
    for folder in expected_folders:
        folder_path = Path(folder)
        if folder_path.exists():
            pdb_files = list(folder_path.glob("*.pdb"))
            print(f"✅ {folder}: {len(pdb_files)} PDB files")
            total_files += len(pdb_files)
        else:
            print(f"❌ {folder}: Missing")
            return False
    
    print(f"✅ Total PDB files: {total_files}")
    return total_files > 0

def check_application():
    """Check if the main application can be imported."""
    print("\n🔧 Checking application...")
    
    try:
        from qsar_web_app import get_ligand_list
        print("✅ Application imports successfully")
        
        # Test ligand detection
        ligands = get_ligand_list('Alpha_CE_Combined')
        if ligands:
            print(f"✅ Ligand detection working: {len(ligands)} ligands found")
            return True
        else:
            print("❌ No ligands detected")
            return False
            
    except Exception as e:
        print(f"❌ Application error: {e}")
        return False

def check_launch_scripts():
    """Check if launch scripts are present."""
    print("\n🚀 Checking launch scripts...")
    
    scripts = ['run_app.bat', 'run_app.ps1']
    for script in scripts:
        if Path(script).exists():
            print(f"✅ {script} - Present")
        else:
            print(f"❌ {script} - Missing")
    
    return True

def main():
    """Run all verification checks."""
    print("=" * 60)
    print("QSAR GUI Setup Verification")
    print("=" * 60)
    
    checks = [
        check_python_version,
        check_dependencies,
        check_data_files,
        check_application,
        check_launch_scripts
    ]
    
    passed = 0
    total = len(checks)
    
    for check in checks:
        if check():
            passed += 1
        print()
    
    print("=" * 60)
    print(f"Verification Results: {passed}/{total} checks passed")
    
    if passed == total:
        print("🎉 All checks passed! Your QSAR GUI is ready to use.")
        print("\nTo start the application:")
        print("  Windows: Double-click run_app.bat")
        print("  PowerShell: .\\run_app.ps1")
        print("  Manual: streamlit run qsar_web_app.py")
        return True
    else:
        print("⚠️  Some checks failed. Please fix the issues above.")
        print("\nCommon solutions:")
        print("  1. Install Python 3.9+: https://python.org")
        print("  2. Install dependencies: pip install -r requirements.txt")
        print("  3. Ensure all data folders are present")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 