ms_table = Arrow.Table("/Users/nathanwamsley/Data/Mar_2025/Kevin_DE_Tag_Pioneer/DE_benchmark_016.arrow")

ms1_scan_idxs = [i for i in range(1, length(ms_table[:msOrder])) if ms_table[:msOrder][i] == 1]
_tasks_, _ = partitionScansToThreadsMS1(
    ms_table[:mz_array],
    ms_table[:msOrder],
    7
)


scans_captured = sort(vcat([task[end] for task in _tasks_]...))
scans_captured = scans_captured[scans_captured.>0]
Set(ms1_scan_idxs)==Set(scans_captured)