$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
$runtimeCommand = Get-Command node -ErrorAction SilentlyContinue
$runtimePath = if ($runtimeCommand) { $runtimeCommand.Source } else { Join-Path $env:USERPROFILE '.cache\codex-runtimes\codex-primary-runtime\dependencies\node\bin\node.exe' }
if (-not (Test-Path -LiteralPath $runtimePath)) {
    Write-Host 'Node.js 20 or newer is required. Install it from https://nodejs.org and run this launcher again.'
    Read-Host 'Press Enter to close'
    exit 1
}
$siteAddress = 'http://127.0.0.1:4173'
try { $runningSite = Invoke-WebRequest -Uri $siteAddress -TimeoutSec 2 -UseBasicParsing } catch { $runningSite = $null }
if ($runningSite -and $runningSite.Content -match 'Polymer Programming') {
    Start-Process $siteAddress
    exit
}
Write-Host "Starting Polymer Programming at $siteAddress"
if ($runningSite) {
    Write-Host 'Port 4173 is already being used by another application. Run server.mjs with a different PORT.'
    Read-Host 'Press Enter to close'
    exit 1
}
$logDirectory = Join-Path $PSScriptRoot '.preview'
New-Item -ItemType Directory -Path $logDirectory -Force | Out-Null
$siteProcess = Start-Process -FilePath $runtimePath -ArgumentList 'server.mjs' -WorkingDirectory $PSScriptRoot -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $logDirectory 'server.log') -RedirectStandardError (Join-Path $logDirectory 'server.error.log')
try {
    $siteReady = $false
    for ($attempt = 0; $attempt -lt 30; $attempt++) {
        if ($siteProcess.HasExited) { throw 'The server could not start. See .preview/server.error.log for details.' }
        try {
            $siteResponse = Invoke-WebRequest -Uri $siteAddress -TimeoutSec 1 -UseBasicParsing
            if ($siteResponse.StatusCode -eq 200) { $siteReady = $true; break }
        } catch { Start-Sleep -Milliseconds 150 }
    }
    if (-not $siteReady) { throw 'The local server did not become ready.' }
    Start-Process $siteAddress
    Write-Host 'Your website is ready. Keep this window open while using it.'
    Read-Host 'Press Enter to stop the server and close this window'
} finally {
    if (-not $siteProcess.HasExited) { Stop-Process -Id $siteProcess.Id -ErrorAction SilentlyContinue }
}
