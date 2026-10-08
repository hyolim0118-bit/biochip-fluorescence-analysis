@echo off
setlocal
cd /d "%~dp0"
if exist ".venv\Scripts\python.exe" (
    ".venv\Scripts\python.exe" "launch_analysis.py"
) else (
    python "launch_analysis.py"
)
pause
