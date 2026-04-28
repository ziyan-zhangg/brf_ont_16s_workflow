# --- Determine next version tag ---
$LAST_TAG = git tag --list "v*" | Sort-Object { [Version]($_.TrimStart('v')) } | Select-Object -Last 1

if ($LAST_TAG) {
    $VER_NUM = [Version]($LAST_TAG.TrimStart("v"))
    $NEW_VER = "v{0}.{1:00}" -f $VER_NUM.Major, ($VER_NUM.Minor + 3)
} else {
    $NEW_VER = "v1.00"
}

$VER = $NEW_VER
$DATE = "Date: {0}" -f (Get-Date -Format "yyyy-MM-dd")

$COMMENT = @"
* Skip some version numbers to make it more clear when a new version is released.
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