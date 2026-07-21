# run_experiments.ps1
# Runs all paper experiments and archives results under a timestamped folder.
# Each invocation creates a fresh Experiments_DDMMYY_HHmm directory — no skipping.

$root     = "C:\Users\heydatey\source\repos\gismo"
$binary   = "$root\build\bin\Release\poissonTHB_example.exe"
$joystick = "$root\filedata\generatedMPs\joystick_approximation_fine_L3_NLO.xml"
$mask     = "$root\filedata\generatedMPs\mask_L3_PatchOnly689_noR.xml"
$hexagon  = "$root\filedata\generatedMPs\hexagon_wiggly.xml"

# ---------------------------------------------------------------------------
# Timestamped root folder for this entire run set
# ---------------------------------------------------------------------------
$timestamp = Get-Date -Format "ddMMyy_HHmm"
$runSetDir = "$root\Dokumentierung\runs\Experiments_$timestamp"
New-Item -ItemType Directory -Path $runSetDir -Force | Out-Null
Write-Host "Run set: $runSetDir"
Write-Host ""

# ---------------------------------------------------------------------------
# Experiment table
# Each entry: name, geometry path (for existence check), extra args
# ---------------------------------------------------------------------------
$experiments = @(
    # Class 1 — tolerance variation (hexagon_wiggly: boundary wave eps=0.03, k=5, no interior distortion)
    # ε_g = 0.10 is loose for the clean interior; ε_f is calibrated to the wave amplitude (0.03).
    # Tight: ε_f < wave amplitude → boundary cells rejected; Medium: ε_f ≈ amplitude; Loose: ε_f >> amplitude.
    @{ name="class1_tight";       geo=$hexagon;  args=@("--epsilon-g","0.10","--epsilon-f","0.020",$hexagon) },
    @{ name="class1_medium";      geo=$hexagon;  args=@("--epsilon-g","0.10","--epsilon-f","0.030",$hexagon) },
    @{ name="class1_loose_f";     geo=$hexagon;  args=@("--epsilon-g","0.10","--epsilon-f","0.050",$hexagon) },

    # Class 2 — projection stage: LS / LO / NLO
    @{ name="class2_ls_only";     geo=$mask;     args=@("--ls-only",$mask) },
    @{ name="class2_lo_only";     geo=$mask;     args=@("--lo-only",$mask) },
    @{ name="class2_full";        geo=$mask;     args=@($mask) },

    # Class 3 — locality parameter lambda, joystick geometry
    @{ name="class3_joy_global";  geo=$joystick; args=@($joystick) },
    @{ name="class3_joy_lambda0"; geo=$joystick; args=@("--local-fitting","--lambda","0",$joystick) },
    @{ name="class3_joy_lambda1"; geo=$joystick; args=@("--local-fitting","--lambda","1",$joystick) },
    @{ name="class3_joy_lambda2"; geo=$joystick; args=@("--local-fitting","--lambda","2",$joystick) },
)

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
if (-not (Test-Path $binary)) {
    Write-Host "ERROR: binary not found: $binary"
    exit 1
}

# ---------------------------------------------------------------------------
# Run loop
# ---------------------------------------------------------------------------
$summary = @()

foreach ($exp in $experiments) {
    $name   = $exp.name
    $outDir = "$runSetDir\$name"

    # Check geometry prerequisite
    if (-not (Test-Path $exp.geo)) {
        $msg = "SKIP $name — geometry missing: $($exp.geo)"
        Write-Host $msg
        $summary += $msg
        continue
    }

    New-Item -ItemType Directory -Path $outDir -Force | Out-Null

    Write-Host "=== $name ===" -ForegroundColor Cyan
    $t0 = Get-Date
    & $binary @($exp.args) *> "$outDir\cmd_output.txt"
    $elapsed = [math]::Round(((Get-Date) - $t0).TotalSeconds, 1)

    # Copy log files the binary writes to the working directory
    Get-ChildItem -Path $root -MaxDepth 1 -Filter "*_summary_log.txt" |
        Copy-Item -Destination $outDir
    Get-ChildItem -Path $root -MaxDepth 1 -Filter "*_logFile_poissonTHB_example.txt" |
        Copy-Item -Destination $outDir

    $msg = "DONE $name — $elapsed s"
    Add-Content "$outDir\timing.txt" $msg
    Write-Host "$msg`n"
    $summary += $msg
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
Write-Host "=== SUMMARY ===" -ForegroundColor Green
$summary | ForEach-Object { Write-Host $_ }
$summary | Set-Content "$runSetDir\summary.txt"
