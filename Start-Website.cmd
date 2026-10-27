@echo off
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0Start-Website.ps1"
if errorlevel 1 pause
