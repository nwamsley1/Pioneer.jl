param(
    [Parameter(Mandatory = $true)]
    [string]$InstallKey
)

$ErrorActionPreference = 'Stop'

function New-PioneerGuid {
    param([Parameter(Mandatory = $true)][string]$Purpose)

    $inputBytes = [System.Text.Encoding]::UTF8.GetBytes("pioneer-installer-v1`:$Purpose`:$InstallKey")
    $digest = [System.Security.Cryptography.SHA256]::HashData($inputBytes)
    $guidBytes = [byte[]]::new(16)
    [Array]::Copy($digest, $guidBytes, 16)

    # Mark the deterministic identifier as an RFC 4122 version-5-style UUID.
    # The bytes come from SHA-256 rather than SHA-1, but the version/variant
    # bits make the generated values well-formed GUIDs for WiX/MSI.
    $guidBytes[7] = ($guidBytes[7] -band 0x0f) -bor 0x50
    $guidBytes[8] = ($guidBytes[8] -band 0x3f) -bor 0x80
    return ([Guid]::new($guidBytes)).ToString().ToUpperInvariant()
}

$values = [ordered]@{
    msi_upgrade_code = New-PioneerGuid 'msi-upgrade'
    bundle_upgrade_code = New-PioneerGuid 'bundle-upgrade'
    permissions_component_guid = New-PioneerGuid 'component-permissions'
    shortcut_component_guid = New-PioneerGuid 'component-shortcut'
    desktop_shortcut_component_guid = New-PioneerGuid 'component-desktop-shortcut'
    path_component_guid = New-PioneerGuid 'component-path'
}

foreach ($entry in $values.GetEnumerator()) {
    "$($entry.Key)=$($entry.Value)" | Out-File -FilePath $env:GITHUB_OUTPUT -Encoding utf8 -Append
    Write-Host "$($entry.Key)=$($entry.Value)"
}
