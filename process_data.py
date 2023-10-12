

from benchmark_list import all_benchmarks
from data_processing_functions import *
from make_figures import *
from pathlib import Path
import re

cvc5_dict = {}

twenty_minutes = 1200

# Sequential data

seq_data = get_sequential_data(benchmarks=all_benchmarks, 
                               file_path="data/runs/sequential_cvc5.csv", 
                               timeout=twenty_minutes)

cvc5_dict[("sequential", "none", "na", 1)] = seq_data

print((f"sequential solved "
       f"{get_num_solved(seq_data)} of "
       f"{len(all_benchmarks)} with par-2 = "
       f"{sum(seq_data)}"))


# Scramble data

scram_log_file = "data/partitioning_or_scrambling/combined_scram_times.csv"
scram_loc = "data/runs/scramble/"
scram_data, scram_log_data = \
    build_scamble_dictionary(benchmarks=all_benchmarks, 
                               scram_logfile=scram_log_file,
                               scram_loc=scram_loc)

current_scram_data = []
for i in range(2, 257):
    current_scram_data = get_scrambled_data(benchmarks=all_benchmarks,
                                    n=i, 
                                    scramdata=scram_data, 
                                    scram_logdata=scram_log_data,
                                    timeout=twenty_minutes,
                                    prev_results=current_scram_data)
    print((f"scram data {i} solved {get_num_solved(current_scram_data)} of "
           f"{len(all_benchmarks)} with par-2 = {sum(current_scram_data)}"))
    cvc5_dict[("scramble", "none", "na", i)] = current_scram_data





# Get all data for a given strategy and set of benchmarks
def get_all_data(s, p, inf, benches, p_type):
    partitioning_times = []
    base_times = []
    total_times = []
    
    info = f"{inf}_{s}_{p}p"
    
    # Need to check the logs to see if, in rare cases, this timed out while 
    # partitioning.
    log_colms=['benchmark', 'solved', 'timed_out', 'result', 'time_real'] 
    log_file = f"data/logs/{info}_logs.csv"
    # print(log_file)
    log_data = pd.read_csv(log_file, names=log_colms)

    # Get the partitioning information for this strategy.
    part_colms =['benchmark', 'num_partitions', 'partitioning_time'] 
    part_file = (f"data/partitioning_or_scrambling/{info}_partitioning_times.csv")
    # print(part_file)
    p_data = pd.read_csv(part_file, names=part_colms)

    for b in benches:
        name = Path(b).stem
        b_log = log_data
        
        if not b_log.empty:
            # If this benchmark is in the log file, then we failed to 
            # partition it. Note that this situation occurs for only a few
            # problems, and none of the benchmarks fail to partition for 
            # smaller numbers of partitions with any strategy. 
            filtered = list(b_log[b_log["benchmark"].str.contains(b)].result)
            if len(filtered) > 0:
                
                # We count failure to partition as a failure to solve.
                rpt = 0
                rbt = 2*twenty_minutes
                rtt = 2*twenty_minutes

                partitioning_times.append(rpt)
                base_times.append(rbt)
                total_times.append(rtt)

                continue 
                
        
        part_time = list(p_data[p_data["benchmark"].str.contains(b)].partitioning_time)
        # print(part_time)
        t = 0

        if len(part_time) > 1:
            sys.exit(f"more than one entry in partitioning time for {b}")
        elif len(part_time) > 0:
            t = part_time[0]
        else:
            sys.exit(f"no data for partitioning time on {b} with strategy {info}")
        partitioning_times.append(t)
        
        
        run_file = (f"data/runs/{info}/{name}_{info}_part.csv")
        run_data = pd.read_csv(run_file)
        
        tts = get_unknown_time(run_data)
        base_times.append(tts)

        # if time to partition plus solve > timeout, we fail.
        if tts == 2*twenty_minutes or (tts + t > twenty_minutes): 
            total_times.append(2*twenty_minutes)
        else:
            total_times.append(tts + t)
        
    return partitioning_times, base_times, total_times


# Partitioning data
ss = ["lemma", "decision", "heap"]
np = [2,4,8,16,32,64,128]

# Cubing data
for s in ss:
    for n in np:
        pt, bt, tt = get_all_data(s, n, "cubing", all_benchmarks, "cubes")
        print("cubing", s, n, "solved", get_num_solved(tt), "of", len(tt), "par-2", sum(tt))
        cvc5_dict[("cubes", s, "na", n)] = tt
print("Cubing data done")

# Climit DNC data
for s in ss:
    for n in np:
        pt, bt, tt = get_all_data(s, n, "climit_dncs", all_benchmarks, "dncs")
        print("climit_dncs", s, n, "solved", get_num_solved(tt), "of", len(tt), "par-2", sum(tt))
        cvc5_dict[("dncs", s, "climit", n)] = tt

print("climit dncs data done")

        
# Time DNC data
for s in ss:
    for n in np:
        pt, bt, tt = get_all_data(s, n, "time_dncs", all_benchmarks, "dncs")
        print("time_dncs", s, n, "solved", get_num_solved(tt), "of", len(tt), "par-2", sum(tt))
        cvc5_dict[("dncs", s, "tlimit", n)] = tt

print("tlimit dncs data done")


###### random cubes
for n in np:
    pt, bt, tt = get_all_data("heap", n, "randcubes", all_benchmarks, "cubes")
    print("randcubes", n, "solved", get_num_solved(tt), "of", len(tt), "par-2", sum(tt))
    cvc5_dict[("cubes", "heap", "RAND", n)] = tt



##### Get all data for a given OpenSMT strategy and set of benchmarks
def get_all_data_osmt(p, inf, benches):
    partitioning_times = []
    base_times = []
    total_times = []
    
    info = f"{inf}_{p}p"
    
    log_file =  f"data/logs/{info}_logs.csv"
    part_file = f"data/partitioning_or_scrambling/{info}_partition_times.csv"
    
    # Need to check the logs to see if solved while partitioning
    log_colms=['benchmark', 'solved', 'timed_out', 'result', 'time_real'] 
    part_colms =['benchmark', 'num_partitions', 'partitioning_time'] 
    
    log_data = pd.read_csv(log_file, names=log_colms)
    part_data = pd.read_csv(part_file, names=part_colms)

    for b in benches:
        name = Path(b).stem

        # OpenSMT2 was allowed up to 21 minutes to partition, so if there is no
        # partitioning data, then we can safely call it a failure. 
        # 21 minutes is the sum of twenty minutes plus one minute to partition,
        # which is what was used to run the cvc5 versions. However, any problem
        # where the sum of the time to partition and the time to solve exceeds
        # twenty minutes is counted as a timeout for both partitioning solvers.
        # Additionally, all benchmarks that timed out during partitioning for 
        # number of partitions n = 2^k also timed out for n = 2^(k+1),
        # so in the interest of time, some problems were skipped for larger
        # values of n for the OpenSMT2 lookahead and deep strategies, which 
        # timed out on dozens of problems, even at smaller values of n. 
        # OpenSMT2's scatter strategy successfully partitioned all 214 
        # problems for all values of n. 
        part_time = list(part_data[part_data["benchmark"].str.contains(b)].partitioning_time)
        t = 0
        if len(part_time) > 0:
            t = part_time[0]
        else:
            partitioning_times.append(0)
            base_times.append(2*twenty_minutes)
            total_times.append(2*twenty_minutes)
            continue

        partitioning_times.append(t)
        
        run_file = (f"data/runs/{info}/{name}_{info}_part.csv")
        
        run_data = pd.read_csv(run_file)
        
        tts = get_unknown_time(run_data)
        
        base_times.append(tts)
        if tts == 2*twenty_minutes or (tts + t > twenty_minutes): 
            total_times.append(2*twenty_minutes)
        else:
            total_times.append(tts + t)
        
    return partitioning_times, base_times, total_times



osmt_dict = {}
osmt_friendly = {}

# Scatter
for n in [2,4,8,16,32,64,128]:
    pt, bt, tt = get_all_data_osmt(n, "osmt_scatter", all_benchmarks)
    osmt_dict[("scatter", n)] = tt
    osmt_friendly[("scatter", "osmt", "osmt", n)] = tt
    print("osmt scatter", n, "solved", get_num_solved(tt), "of", len(tt), "par-2", sum(tt))
    
# Lookahead
for n in [2,4,8,16,32,64,128]:
    pt, bt, tt = get_all_data_osmt(n, "osmt_lookahead", all_benchmarks)
    osmt_dict[("lookahead", n)] = tt
    osmt_friendly[("lookahead", "osmt", "osmt", n)] = tt
    print("osmt lookahead", n, "solved", get_num_solved(tt), "of", len(tt), "par-2", sum(tt))
    
# Lookahead deep
for n in [2,4,8,16,32,64,128]:
    pt, bt, tt = get_all_data_osmt(n, "osmt_deep", all_benchmarks)
    osmt_dict[("deep", n)] = tt
    osmt_friendly[("deep", "osmt", "osmt", n)] = tt
    print("osmt deep", n, "solved", get_num_solved(tt), "of", len(tt), "par-2", sum(tt))
    
make_figure_one(cvc5_dict)
make_figure_two(cvc5_dict)
make_figure_three(cvc5_dict, osmt_dict)
make_figure_four(cvc5_dict, osmt_dict, osmt_friendly)

def get_parallel_runtime(data, cores):
    tasks = list(zip(data['benchmark'], data['time_real'], data['version'], data["result"]))

    core_queue = PriorityQueue()
    for i in range(cores):
        core_queue.put((0, 0, ""))

    finish_times = {}
    sat_finish_times = {}
    result_map = {}
    for v in list(data['version']):
        result_map[v] = "solved"
        

    last_core = None
    while tasks:
        finish_time, core, core_version = core_queue.get()

        if core_version != "":
            if core_version in finish_times.keys():
                finish_times[core_version] = max(finish_times[core_version], finish_time)
            else:
                finish_times[core_version] = finish_time

        benchmark, time_real, version, result = tasks.pop(0)
        
        finish_time = finish_time + time_real
        if version != "":
            if version not in finish_times.keys():
                if result == "sat":
                    sat_finish_times[version] = finish_time
                    finish_times[version] = finish_time
                    result_map["version"] = "solved"
                elif result == "unknown":
                    result_map["version"] = "unsolved"
                    finish_times[version] = finish_time
                else:
                    finish_times[version] = finish_time
            elif version in finish_times.keys():
                if result == "sat":
                    finish_times[version] = min(finish_times[version], finish_time)
                    if version in sat_finish_times.keys():
                        sat_finish_times[version] = min(sat_finish_times[version], finish_time)
                    else:
                        sat_finish_times[version] = finish_time
                    result_map["version"] = "solved"
                elif result == "unknown":
                    result_map["version"] = "unsolved"
                    finish_times[version] = max(finish_times[version], finish_time) 
                else:
                    finish_times[version] = max(finish_times[version], finish_time) 
                    
        core_queue.put((finish_time, core, version))
        last_core = (finish_time, core, version)
    
    finish_time, core, version = last_core
    
    if version != "":
        if version not in finish_times.keys():
            if result == "sat":
                sat_finish_times[version] = finish_time
                finish_times[version] = finish_time
                result_map["version"] = "solved"
            elif result == "unknown":
                result_map["version"] = "unsolved"
                finish_times[version] = finish_time
            else:
                finish_times[version] = finish_time
        elif version in finish_times.keys():
            if result == "sat":
                finish_times[version] = min(finish_times[version], finish_time)
                if version in sat_finish_times.keys():
                    sat_finish_times[version] = min(sat_finish_times[version], finish_time)
                else:
                    sat_finish_times[version] = finish_time
                result_map["version"] = "solved"
            elif result == "unknown":
                result_map["version"] = "unsolved"
                finish_times[version] = max(finish_times[version], finish_time) 
            else:
                finish_times[version] = max(finish_times[version], finish_time) 
    
    if "solved" in [result_map[k] for k, v in finish_times.items()]:
        if len(list(sat_finish_times.keys())) == 0:
            total_t = min([v for k, v in finish_times.items() if result_map[k] == "solved"])
            if total_t > 1200:
                return 2400
            else:
                return total_t
        else:
            total_t = min([v for k, v in sat_finish_times.items() if result_map[k] == "solved"])
            if total_t > 1200:
                return 2400
            else:
                return total_t
    else: 
        return 2400


def ser_sort(ser):
    return ser.apply(lambda x: int(re.search(r'_p(\d+)\.smt2', x).group(1)))

def get_all_data_mega(s, p, inf, b, p_type, check_reruns=True):
    partition_times = []
    partition_results = []
    
    info = f"{inf}_{s}_{p}p"
    
    log_colms=['benchmark', 'solved', 'timed_out', 'result', 'time_real'] 
    log_file = f"data/logs/{info}_logs.csv"
    log_data = pd.read_csv(log_file, names=log_colms)
    
    part_colms =['benchmark', 'num_partitions', 'partitioning_time'] 
    part_file = f"data/partitioning_or_scrambling/{info}_partitioning_times.csv"
    p_data = pd.read_csv(part_file, names=part_colms)
    
    name = Path(b).stem
    b_log = log_data

    if not b_log.empty:
        filtered = list(b_log[b_log["benchmark"].str.contains(b)].result)
        if len(filtered) > 0:
            return partition_times, partition_results

    part_time = list(p_data[p_data["benchmark"].str.contains(b)].partitioning_time)
    t = part_time[0]
    
    run_file = (f"data/runs/{info}/{name}_{info}_part.csv")
    
    run_data = pd.read_csv(run_file)
    run_data = run_data.sort_values(by='benchmark', key=ser_sort)

    for p in list(run_data["time_real"]):
        partition_times.append(p + t)
    for r in list(run_data["result"]):
        partition_results.append(r)
        
    return partition_times, partition_results


def get_all_data_mega_osmt(s, p, inf, b):
    partition_times = []
    partition_results = []
    
    info = f"{inf}_{p}p"
    
    log_file =  f"data/logs/{info}_logs.csv"
    part_file = f"data/partitioning_or_scrambling/{info}_partition_times.csv"
    
    log_colms=['benchmark', 'solved', 'timed_out', 'result', 'time_real'] 
    part_colms =['benchmark', 'num_partitions', 'partitioning_time'] 
    
    log_data = pd.read_csv(log_file, names=log_colms)
    part_data = pd.read_csv(part_file, names=part_colms)

    name = Path(b).stem

    b_log = log_data

    if not b_log.empty:
        filtered = list(b_log[b_log["benchmark"].str.contains(b)].result)
        if len(filtered) > 0:
            return partition_times, partition_results

    p_data = part_data

    part_time = list(p_data[p_data["benchmark"].str.contains(b)].partitioning_time)
    t = part_time[0]

    run_file = (f"data/runs/{info}/{name}_{info}_part.csv")
                 
    run_data = pd.read_csv(run_file)
    run_data = run_data.sort_values(by='benchmark', key=ser_sort)
    
    for p in list(run_data["time_real"]):
        partition_times.append(p + t)
    for r in list(run_data["result"]):
        partition_results.append(r)
        
    return  partition_times, partition_results


scheduling_x = range(2,257)
scheduling_y = []

for n in scheduling_x:
    lst = get_recommended_graduated_portfolio(n, cvc5_dict, osmt_friendly)
    tts = []
    i = 0
    for b in all_benchmarks:
        benches = []
        results = []
        times = []
        versions = []
        np = [2, 4, 8, 16, 32, 64, 128]
        for x in np:
            pts, prs = get_all_data_mega_osmt("opensmt_partitioning",  x, "osmt_scatter", b)
            if len(pts) > 0:
                benches = benches + [b]*len(prs)
                results = results + prs
                times = times + pts
                versions = versions + [f"osmt_{x}_scattering_{b}"]*len(prs)

            pts, prs = get_all_data_mega("decision", x, "cubing", b, "cubes", check_reruns=True)
            if len(pts) > 0:
                benches = benches + [b]*len(prs)
                results = results + prs
                times = times + pts
                versions = versions + [f"decision_{x}_cubing_{b}_cubes"]*len(prs)

        data = pd.DataFrame({
            'benchmark': benches,
            'result': results,
            'time_real': times,
            'version': versions
        })

        ti = get_parallel_runtime(data, n)
        tts.append(ti)
        i+=1
    
    scheduling_y.append(sum(tts))


hybrid_scheduling_x = range(4,257)
hybrid_scheduling_y = []

for n in hybrid_scheduling_x:
    tts = []
    i = 0
    for b in all_benchmarks:
        benches = []
        results = []
        times = []
        versions = []
        np = [2, 4, 8, 16, 32, 64, 128]
        for x in np:
            pts, prs = get_all_data_mega_osmt("opensmt_partitioning",  x, "osmt_scatter", b)
            if len(pts) > 0:
                benches = benches + [b]*len(prs)
                results = results + prs
                times = times + pts
                versions = versions + [f"osmt_{x}_scattering_{b}"]*len(prs)

            pts, prs = get_all_data_mega("decision", x, "cubing", b, "cubes", check_reruns=True)
            if len(pts) > 0:
                benches = benches + [b]*len(prs)
                results = results + prs
                times = times + pts
                versions = versions + [f"decision_{x}_cubing_{b}_cubes"]*len(prs)

        data = pd.DataFrame({
            'benchmark': benches,
            'result': results,
            'time_real': times,
            'version': versions
        })

        ti = get_parallel_runtime(data, n//2)
        
        tts.append(ti)
        i+=1
    
    scram_tts = cvc5_dict[("scramble", "none", "na", n-(n//2))]
    hybrid_scheduling_y.append(sum(get_total_time([tts, scram_tts])))


np = [2,4,8,16,32,64,128,256]
plt.scatter([1], [sum(cvc5_dict[("sequential", "none", "na", 1)])], color="black", label="sequential",
           marker="X")

scramY = [sum(cvc5_dict[("scramble", "none", "na", i)]) for i in range(2,257)]
plt.plot(range(2,257), scramY, label="scrambling portfolio", color="black", linestyle=":", alpha=.5)

PPX = [get_req_num_cores(get_recommended_graduated_portfolio(n, cvc5_dict, osmt_friendly)) for n in np]
PPY = []
for n in np:
    lst = get_recommended_graduated_portfolio(n, cvc5_dict, osmt_friendly)
    tt = []
    for l in lst:
        if l in cvc5_dict.keys():
            tt.append(cvc5_dict[l])
        else:
            tt.append(osmt_friendly[l])
    print("RP", n, get_num_solved(get_total_time(tt)))
    PPY.append(sum(get_total_time(tt)))

plt.plot(PPX, PPY, label="recommended portfolio", color="black", 
          marker="^", markeredgecolor="black", 
            markerfacecolor="white", alpha=.5, markersize=7)
rspX = [get_req_num_cores(get_rec_scrambling(n, cvc5_dict, osmt_friendly)) for n in np]
rspY = []
for n in np:
    lst = get_rec_scrambling(n, cvc5_dict, osmt_friendly)
    tt = []
    for l in lst:
        if l in cvc5_dict.keys():
            tt.append(cvc5_dict[l])
        else:
            tt.append(osmt_friendly[l])
    print("hybrid", n, get_num_solved(get_total_time(tt)), sum(get_total_time(tt)))
    rspY.append(sum(get_total_time(tt)))

plt.plot(rspX, rspY, label="hybrid",
         alpha=.5, color="black", marker="o")

plt.plot(scheduling_x, scheduling_y, label="recommended portfolio (multijob)", linestyle="--",
         color="black", alpha=.5)

plt.plot(hybrid_scheduling_x, hybrid_scheduling_y, label="hybrid (multijob)", 
         color="black")

plt.legend()
plt.xscale("function", functions=(forward, inverse))
plt.xticks([1,2,4,8,16,32,64,128,256])

plt.xlabel("Number of Cores")
plt.ylabel("PAR-2 Score")
plt.savefig("figures/figure_5.png", bbox_inches='tight')

