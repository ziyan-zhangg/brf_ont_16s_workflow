# --- Determine next version tag ---
$LAST_TAG = git tag --list "v*" | Sort-Object { [Version]($_.TrimStart('v')) } | Select-Object -Last 1

if ($LAST_TAG) {
    $VER_NUM = [Version]($LAST_TAG.TrimStart("v"))
    $NEW_VER = "v{0}.{1:00}" -f $VER_NUM.Major, ($VER_NUM.Minor + 1)
} else {
    $NEW_VER = "v1.00"
}

$VER = $NEW_VER
$DATE = "Date: {0}" -f (Get-Date -Format "yyyy-MM-dd")

$COMMENT = @"
* Fix: now prep script will look into fastq_pass and fastq_fail folders for fastq files, instead of looking into the root folder. 
  This is to accomodate the strict filtering applied for plasmid.
"@

Write-Host "Using version: $VER"

# --- Update changelog ---
Move-Item -Path "changelog.txt" -Destination "changelog_old.txt"
Set-Content -Path "changelog.txt" -Value $VER
Add-Content -Path "changelog.txt" -Value $DATE
Add-Content -Path "changelog.txt" -Value $COMMENT
Add-Content -Path "changelog.txt" -Value ""
Get-Content -Path "changelog_old.txt" | Add-Content -Path "changelog.txt"
Remove-Item -Path "changelog_old.txt"

# --- Git operations ---
git add .
git commit -m $COMMENT
git tag -a $VER -m $COMMENT

git push -u origin main
git push origin $VER