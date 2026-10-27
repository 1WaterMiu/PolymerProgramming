$ErrorActionPreference = 'Stop'
$pidFile = Join-Path $PSScriptRoot '.preview/server.pid'
if (-not (Test-Path -LiteralPath $pidFile)) { Write-Host 'No managed server is running.'; exit }
$serverProcessId = [int](Get-Content -LiteralPath $pidFile)
$serverProcess = Get-CimInstance Win32_Process -Filter "ProcessId = $serverProcessId"
$serverPath = Join-Path $PSScriptRoot 'server.mjs'
if ($serverProcess -and $serverProcess.Name -eq 'node.exe' -and $serverProcess.CommandLine.Contains($serverPath)) {
    Stop-Process -Id $serverProcessId
    Write-Host 'Polymer Programming has stopped.'
} elseif ($serverProcess) {
    throw 'The saved process is not this website server. It has not been stopped.'
}
Remove-Item -LiteralPath $pidFile
