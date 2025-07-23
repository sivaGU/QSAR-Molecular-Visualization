Write-Host "Starting QSAR Molecular Visualization Tool..." -ForegroundColor Green
Write-Host ""
Write-Host "Installing dependencies..." -ForegroundColor Yellow
pip install -r requirements.txt
Write-Host ""
Write-Host "Launching application..." -ForegroundColor Green
streamlit run qsar_web_app.py 