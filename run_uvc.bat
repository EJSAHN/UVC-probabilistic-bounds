@echo off
setlocal EnableExtensions
chcp 65001 >nul

rem =========================================================
rem run_uvc.bat  (GitHub-safe runner)
rem Usage:
rem   run_uvc.bat  "C:\path\to\DATA_DIR"   [OUT_DIR]   [ENV_NAME]
rem
rem Examples:
rem   run_uvc.bat "C:\projects\UV_data"
rem   run_uvc.bat "C:\projects\UV_data" "C:\projects\UV_2\outputs" sap-kmer-gwas
rem
rem Notes:
rem - DATA_DIR must contain Coffee.xlsx and Cacao.xlsx
rem - OUT_DIR default: .\outputs
rem - ENV_NAME optional; if provided we call 'conda activate ENV_NAME'
rem =========================================================

if "%~1"=="" (
  echo [ERROR] Missing DATA_DIR argument.
  echo Usage: run_uvc.bat "C:\path\to\DATA_DIR" [OUT_DIR] [ENV_NAME]
  exit /b 1
)

set "DATA_DIR=%~1"
set "OUT_DIR=%~2"
set "ENV_NAME=%~3"

if "%OUT_DIR%"=="" set "OUT_DIR=%CD%\outputs"

if not exist "%DATA_DIR%" (
  echo [ERROR] DATA_DIR not found: "%DATA_DIR%"
  exit /b 1
)

if not "%ENV_NAME%"=="" (
  call conda activate "%ENV_NAME%"
  if errorlevel 1 (
    echo [ERROR] conda activate failed: "%ENV_NAME%"
    exit /b 1
  )
)

echo [INFO] DATA_DIR = "%DATA_DIR%"
echo [INFO] OUT_DIR  = "%OUT_DIR%"

python uvc_run_all.py --data-dir "%DATA_DIR%" --out-dir "%OUT_DIR%"
if errorlevel 1 (
  echo [ERROR] uvc_run_all.py failed.
  exit /b 1
)

echo [DONE] Outputs written to "%OUT_DIR%"
endlocal
exit /b 0
