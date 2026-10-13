using Test

@testset "storage kind detection" begin
    # Live probe on whatever this machine is: must return a known symbol and never throw.
    @test Pioneer.storage_kind(@__DIR__) in Pioneer.STORAGE_KINDS
    @test Pioneer.storage_kind(joinpath(@__DIR__, "does_not_exist_$(rand(UInt32))")) === :unknown
    @test !Pioneer.storage_is_slow(:ssd) && !Pioneer.storage_is_slow(:network) && !Pioneer.storage_is_slow(:unknown)
    @test Pioneer.storage_is_slow(:hdd) && Pioneer.storage_is_slow(:optical)

    @testset "macOS diskutil plist" begin
        @test Pioneer._macos_plist_kind("<key>SolidState</key>\n\t<true/>") === :ssd
        @test Pioneer._macos_plist_kind("<key>SolidState</key>\n\t<false/>") === :hdd
        @test Pioneer._macos_plist_kind("<key>OpticalMediaType</key><string>DVD-ROM</string>") === :optical
        @test Pioneer._macos_plist_kind("<key>Internal</key><true/>") === :unknown
    end

    @testset "Windows probe output" begin
        @test Pioneer._windows_probe_kind("HDD\r\n") === :hdd
        @test Pioneer._windows_probe_kind("SSD") === :ssd
        @test Pioneer._windows_probe_kind("SCM") === :ssd
        @test Pioneer._windows_probe_kind("OPTICAL") === :optical
        @test Pioneer._windows_probe_kind("NETWORK") === :network
        @test Pioneer._windows_probe_kind("Unspecified") === :unknown
        @test occursin("DriveLetter D", Pioneer._windows_probe_script('D'))
    end

    @testset "Linux mountinfo" begin
        mi = """
        22 1 8:2 / / rw,relatime shared:1 - ext4 /dev/sda2 rw
        30 22 8:17 / /data rw,relatime shared:2 - ext4 /dev/sdb1 rw
        31 22 0:45 / /mnt/ris rw,relatime shared:3 - cifs //storage1/share rw
        32 22 8:33 / /mnt/my\\040disk rw,relatime shared:4 - xfs /dev/sdc1 rw
        """
        @test Pioneer._linux_mount_for("/data/run1.arrow", mi) == ("ext4", "/dev/sdb1")
        @test Pioneer._linux_mount_for("/database/x", mi) == ("ext4", "/dev/sda2")   # /data is not a prefix of /database
        @test Pioneer._linux_mount_for("/mnt/ris/a.arrow", mi) == ("cifs", "//storage1/share")
        @test Pioneer._linux_mount_for("/mnt/my disk/a", mi) == ("xfs", "/dev/sdc1")
        @test Pioneer._linux_major(UInt64(0x0811)) == 8 && Pioneer._linux_minor(UInt64(0x0811)) == 17
    end

    @testset "one probe per volume, links resolved" begin
        mktempdir() do dir
            f1, f2 = joinpath(dir, "a.arrow"), joinpath(dir, "b.arrow")
            write(f1, "x"); write(f2, "y")
            link = joinpath(dir, "link.arrow"); symlink(f1, link)
            @test Pioneer._storage_volume_key(f1) == Pioneer._storage_volume_key(f2) == Pioneer._storage_volume_key(link)
            # a volume of this machine is never reported twice
            @test length(Pioneer.slow_storage_paths([f1, f2, link])) <= 1
        end
    end
end
