# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# Byte layout constants. Offsets are 0-based.

# --- File preamble / device directory ---
const PREAMBLE_MARKER       = 0xa101   # u16 at file offset 0
const VERSION_OFFSET        = 0x24     # u32 format version: 62, 63, 64 or 66
# The directory position varies with the length of the preceding text fields. It is found as the unique position k
# where the u32 pair at k - count_back (version-dependent, below) is (n, n) and the n 16-byte entries at k
# (u32 kind, u32 index, u64 offset) are well formed.
const DEVICE_DIR_ENTRY_SIZE = 16
const DEVICE_DIR_MAX        = 32
const HEAD_SEARCH_BYTES     = 65536
const DEVICE_BLOCK_SIGNATURE = (0x00000001, 0x00000000, 0x00000001)
const KIND_MASS_SPEC        = 0x00000000

# --- Device block ---
const DEV_FIRST_SCAN   = 0x08          # u32
const DEV_LAST_SCAN    = 0x0c          # u32
const DEV_FIRST_RT     = 0x48          # f64 minutes, == RT of first scan
const DEV_LAST_RT      = 0x50          # f64 minutes, == RT of last scan

"""
Everything that differs between format versions: record sizes and the positions of offsets and fields.
"""
struct VersionLayout
    count_back::Int        # directory: distance from the (n, n) count pair to the first entry
    dev_offset_width::Int  # device-block offsets: 8 = u64, 4 = u32
    dev_scan_table::Int    # device-block positions of: scan table, packet start, region C, region D, self
    dev_packet_start::Int
    dev_region_c::Int
    dev_region_d::Int
    dev_self::Int
    scan_record::Int       # bytes per scan-table record
    st_pkt_off::Int        # packet offset position in the record (relative to packet start)
    st_pkt_off_width::Int  # 8 = u64, 4 = u32
    st_cycle::Int          # instrument cycle counter position, -1 if absent
    rc_fixed::Int          # region C: fixed part, precursor entry, tail
    rc_prec::Int
    rc_tail::Int
    peak_words::Int        # 4-byte words per peak in the peak list: 2 = (f32 m/z, f32 int); 3 = (f64 m/z, f32 int)
    dev_region_a::Int      # device-block positions of region A (logged records) and region B (typed field lists)
    dev_region_b::Int
end
# version 66
const LAYOUT_V66 = VersionLayout(0x314, 8, 0x1cf0, 0x1cf8, 0x1d18, 0x1d20, 0x1d30, 88, 72, 8, 80, 136, 56, 8, 2, 0x1d00, 0x1d08)
# version 64
const LAYOUT_V64 = VersionLayout(0x314, 8, 0x1cf0, 0x1cf8, 0x1d18, 0x1d20, 0x1d30, 80, 72, 4, 76, 128, 32, 4, 2, 0x1d00, 0x1d08)
# version 63; version 62 is as 63 but with a 120-byte region C fixed part and 12-byte peaks (f64 m/z, f32 intensity)
const LAYOUT_V63 = VersionLayout(8, 4, 0x1c, 0x20, 0x1cc8, 0x1ccc, 0x1ce4, 72, 0, 4, -1, 128, 32, 4, 2, 0x24, 0x28)
const LAYOUT_V62 = VersionLayout(8, 4, 0x1c, 0x20, 0x1cc8, 0x1ccc, 0x1ce4, 72, 0, 4, -1, 120, 32, 4, 3, 0x24, 0x28)
const LAYOUTS = Dict(66 => LAYOUT_V66, 64 => LAYOUT_V64, 63 => LAYOUT_V63, 62 => LAYOUT_V62)

# --- Scan table: one record per scan (size and packet-offset position by version, above) ---
const ST_INDEX     = 4     # u32, 0-based
const ST_MSCLASS   = 8     # u32, position in the acquisition cycle: 0 MS1, 1.. MS2 (0xffffffff when not stored)
const ST_SCAN      = 12    # u32, scan number
const ST_PKT_TYPE  = 16    # u32 packet type: 21 profile + centroids, 20 centroids (high-res), 19 profile only (ion trap), 18 centroids (ion trap)
const ST_PKT_SIZE  = 20    # u32, packet size in bytes
const ST_RT        = 24    # f64, minutes
const ST_TIC       = 32    # f64, = sum of peak intensities
const ST_BP_INT    = 40    # f64, intensity of highest peak
const ST_BP_MZ     = 48    # f64, m/z of highest peak
const ST_LOW_MZ    = 56    # f64
const ST_HIGH_MZ   = 64    # f64

# --- Region C: one variable-length record per scan, after a 4-byte lead ---
# record = [136-byte fixed part][u32 n_prec][n_prec x 56][u32 n_ranges][n_ranges x 16][u32 n_coef][n_coef x 8]
#          [u32 n_extra][n_extra x 8][8 bytes]; the last record ends exactly where region D starts.
const RC_LEAD       = 4     # u32: 0 (version 66) or the record count (versions 62-64)
const RC_RANGE_SIZE = 16    # fixed / precursor-entry / tail sizes are in VersionLayout
const RC_MS_ORDER   = 6     # u8 in the fixed part
const RC_ANALYSER   = 40    # u8 analyser code: 0 ion trap, 4 orbitrap, 7 astral analyser
const RC_ENERGY     = 112   # f64 copy of the collision energy (zero in some files)
# precursor entry (56 bytes), offsets within the entry
const PE_CENTER = 0         # f64 isolation centre (MS1: centre of the scan range, when present)
const PE_WIDTH  = 8         # f64 isolation width
const PE_ENERGY = 16        # f64 collision energy
# first range entry: f64 low, f64 high (== scan-table range)

# --- Spectrum packet ---
const PKT_HEADER     = 40  # bytes
const PKT_W_PROFILE  = 4   # u32 word count of the profile section
const PKT_W_PEAKS    = 8   # u32 word count of the peak list (= 1 + peak_words × n)
const PKT_W_REST     = (16, 20, 24, 28)  # u32 word counts of four trailing sections
const PKT_LOW_MZ     = 32  # f32
const PKT_HIGH_MZ    = 36  # f32
# Section 3 (W3 == n words): one u32 per peak
const PW_INDEX_MASK  = 0x0007ffff   # bits 0-18: peak index
const PW_CHARGE_SHIFT = 24          # bits 24-31: charge (0 = none); bits 19-23: flags (21 = reference ion)

# packet types with a decoded layout; anything else is refused.
# 0x13: ion-trap profile-only scans (profile points, no peak list: W2 = 0)
const KNOWN_PACKET_TYPES = (0x00000012, 0x00000013, 0x00000014, 0x00000015)
const ANALYSER_NAMES = Dict(0x00 => "ion trap", 0x04 => "orbitrap", 0x07 => "astral")
