from queue import PriorityQueue
import matplotlib.pyplot as plt
from data_processing_functions import *
import re

def get_total_time(lst):
    total_times = [2400] * len(lst[0])
    for sub_l in lst:
        for i in range(len(sub_l)):
            if total_times[i] > sub_l[i]:
                total_times[i] = sub_l[i]
    return total_times

def get_portfolio1(n):
    if n == 2:
        lst = [("cubes", "decision", "na", n)]
    elif n == 4:
        lst = [("cubes", "decision", "na", n/2), ("cubes", "lemma", "na", n/2)]
    else:
        lst = [("cubes", "decision", "na", n/2),
              ("cubes", "lemma", "na", n/2)]
    return lst

def get_portfolio2(n):
    if n == 2:
        lst = [("cubes", "decision", "na", n)]
    elif n == 4:
        lst = [("cubes", "decision", "na", n/2), ("cubes", "lemma", "na", n/2)]
    else:
        lst = [("cubes", "decision", "na", n/2),
              ("cubes", "lemma", "na", n/4),
              ("scatter", "osmt", "osmt", n/4)]
    return lst

def get_req_num_cores(lst: list):
    if len(lst) == 0:
        return 0
    su = 0
    seen_lst = []
    for strat in lst:
        su += strat[3] # 4th element is number of partitions for cvc5 data
    return su

def get_req_num_cores_osmt(lst: list):
    if len(lst) == 0:
        return 0
    su = 0
    seen_lst = []
    for strat in lst:
        su += strat[1] # 1st element is number of partitions for osmt data 
    return su

def get_scatter_portfolio_osmt(n, osmt_dict):
    list_of_keys = [k for k,v in osmt_dict.items() if k[0] == "scatter"]
    sorted_list = sorted(list_of_keys, key=lambda dat: dat[1])
    sorted_list = sorted(sorted_list, key=lambda dat: dat[0])
    current_n = 0
    lst = []
    i = 0
    if n == 2:
        return [("scatter", n)]
    while current_n < n:
        lst.append(sorted_list[i])
        current_n = get_req_num_cores_osmt(lst)
        i += 1
        if current_n > n:
            lst.pop()
            break
    return lst

# TODO: rename this!!!
def get_graduated_portfolio_v1(n, cvc5_dict, osmt_friendly):
    list_of_keys1 = [k for k,v in osmt_friendly.items() if k[0] == "scatter"]
    list_of_keys2 = [k for k,v in cvc5_dict.items() if (k[0] == "cubes" and k[1] == "decision") or 
                    (k[0] == "cubes" and k[1] == "lemma")]
    list_of_keys = list_of_keys1 + list_of_keys2
    
    sorted_list = sorted(list_of_keys, key=lambda dat: dat[1])
    sorted_list = sorted(sorted_list, key=lambda dat: dat[0], reverse=True)
    sorted_list = sorted(sorted_list, key=lambda dat: dat[3])
    current_n = 0
    lst = []
    i = 0
    if n == 2:
        return [("scatter", "osmt", "osmt", n)]
    while current_n < n:
        lst.append(sorted_list[i])
        current_n = get_req_num_cores(lst)
        i += 1
        if current_n > n:
            lst.pop()
            break
    return lst

# TODO: RENAME THIS!
def get_recommended_graduated_portfolio(n, cvc5_dict, osmt_friendly):
    list_of_keys1 = [k for k,v in osmt_friendly.items() if k[0] == "scatter"]
    list_of_keys2 = [k for k,v in cvc5_dict.items() if (k[0] == "cubes" and k[1] == "decision")]
    list_of_keys = list_of_keys1 + list_of_keys2
    
    sorted_list = sorted(list_of_keys, key=lambda dat: dat[1])
    sorted_list = sorted(sorted_list, key=lambda dat: dat[0], reverse=True)
    sorted_list = sorted(sorted_list, key=lambda dat: dat[3])

    current_n = 0
    lst = []
    i = 0
    if n == 2:
        return [("scatter", "osmt", "osmt", n)]
    while current_n < n:
        lst.append(sorted_list[i])
        current_n = get_req_num_cores(lst)
        i += 1
        if current_n > n:
            lst.pop()
            break
    return lst

def get_rec_scrambling(n, cvc5_dict, osmt_friendly):
    if n == 2:
        lst = get_recommended_graduated_portfolio(n, cvc5_dict, osmt_friendly)
    else:
        lst = get_recommended_graduated_portfolio(n/2, cvc5_dict, osmt_friendly)
        lst.append(("scramble", "none", "na", n-get_req_num_cores(lst)))
    return lst 



def make_figure_one(cvc5_dict):
    ss = ["lemma", "decision", "heap"]
    np = [2, 4, 8, 16, 32, 64, 128]
    
    marker_dict = {
        "lemma": "o",
        "decision": "^",
        "heap": "d"
    }

    label_dict = {
        "lemma": ("check-cl-scatter-spec", "time-cl-scatter-spec"),
        "decision": ("check-decision-scatter-spec", "time-decision-scatter-spec"),
        "heap": ("check-heap-scatter-spec", "time-heap-scatter-spec")
    }
    
    facecolor_dict = {
        "climit": "white",
        "tlimit": "black"
    }

    for s in ss:
        m = marker_dict[s]

        for limit_type, facecolor in facecolor_dict.items():
            xs = []
            ys = []
            for n in np:
                xs.append(n)
                ys.append(sum(cvc5_dict[("dncs", s, limit_type, n)]))
            plt.plot(xs, ys, linestyle="--", color="black", label=label_dict[s][0 if limit_type == "climit" else 1], 
                     marker=m, markeredgecolor="black", markerfacecolor=facecolor)

        plt.legend()
        plt.xlabel("Number of Partitions")
        plt.ylabel("PAR-2 Score")
        plt.xscale("function", functions=(forward, inverse))
        plt.xticks(np)
        plt.savefig(f"figures/figure_1_{s}.png", bbox_inches='tight')
        plt.clf()

def make_figure_two(cvc5_dict):
    np = [2, 4, 8, 16, 32, 64, 128]

    # Plotting for RAND
    randY = [sum(cvc5_dict[("cubes", "heap", "RAND", n)]) for n in np]
    plt.plot(np, randY, linestyle="--", color="black", label=f"time-heap-cube-rand", alpha=.5)

    attribute_dict = {
        "lemma": {
            "marker": "o",
            "label": ("time-cl-cube-spec", "time-cl-scatter-spec")
        },
        "decision": {
            "marker": "s",
            "label": ("time-decision-cube-spec", "time-decision-scatter-spec")
        },
        "heap": {
            "marker": "^",
            "label": ("time-heap-cube-spec", "time-heap-scatter-spec")
        }
    }
    
    msz = 6 ** 2
    data_keys = [("cubes", "na"), ("dncs", "tlimit")]
    face_colors = ["black", "white"]

    for s in attribute_dict.keys():
        for idx, (data_type, limit) in enumerate(data_keys):
            ys = [sum(cvc5_dict[(data_type, s, limit, n)]) for n in np]
            plt.scatter(np, ys, label=attribute_dict[s]["label"][idx], marker=attribute_dict[s]["marker"], 
                        edgecolors="black", facecolors=face_colors[idx], alpha=.5, s=msz)

    plt.xlabel("Number of Partitions")
    plt.ylabel("PAR-2 Score")
    plt.xscale("function", functions=(forward, inverse))
    plt.xticks(np)
    plt.legend()
    plt.savefig("figures/figure_2.png", bbox_inches='tight')
    plt.clf()


def make_figure_three(cvc5_dict, osmt_dict):
    np = [2, 4, 8, 16, 32, 64, 128]
    
    # Define plot configurations in a list of dictionaries
    plot_configs = [
        {"data": osmt_dict, "key": ("scatter",), "label": "OpenSMT2 Scatter", "linestyle": "--", "marker": "X"},
        {"data": osmt_dict, "key": ("lookahead",), "label": "OpenSMT2 Lookahead", "linestyle": "--", "marker": "d"},
        {"data": osmt_dict, "key": ("deep",), "label": "OpenSMT2 alt. Lookahead", "linestyle": "--", "marker": "D"},
        {"data": cvc5_dict, "key": ("cubes", "decision", "na"), "label": "decision-cube", "linestyle": ":", "marker": "s", "facecolor": "white"},
        {"data": cvc5_dict, "key": ("cubes", "lemma", "na"), "label": "cl-cube", "linestyle": ":", "marker": "P", "facecolor": "white"},
        {"data": cvc5_dict, "key": ("dncs", "decision", "tlimit"), "label": "decision-scatter", "linestyle": ":", "marker": "o", "facecolor": "white"}
    ]

    fig, ax = plt.subplots()

    handles = []

    for config in plot_configs:
        ys = [sum(config["data"][config["key"] + (n,)]) for n in np]
        handle, = ax.plot(np, ys, label=config["label"], linestyle=config.get("linestyle", ""),
                          marker=config.get("marker", ""), color="black", alpha=.5, markersize=7,
                          markeredgecolor="black", markerfacecolor=config.get("facecolor", "black"))
        handles.append(handle)

    ax.legend(handles=handles, loc="upper right")
    plt.xlabel("Number of Partitions")
    plt.ylabel("PAR-2 Score")
    plt.xscale("function", functions=(forward, inverse))
    plt.xticks(np)
    plt.savefig("figures/figure_3.png", bbox_inches='tight')
    plt.clf()

def make_figure_four(cvc5_dict, osmt_dict, osmt_friendly):
    np = [2,4,8,16,32,64,128]

    osmt_scatterY = [sum(osmt_dict[("scatter", n)]) for n in np]
    plt.plot(np, osmt_scatterY, label="OpenSMT2 scatter", linestyle="--", marker="X", 
             color="black", alpha=.5)

    p1Y = []
    for n in np:
        my_lst = get_portfolio1(n)
        tt = get_total_time([cvc5_dict[l] for l in my_lst])
        p1Y.append(sum(tt))

    plt.plot(np, p1Y, label="portfolio 1",
             alpha=.5, color="black", marker="d", linestyle="--", markerfacecolor="white", markeredgecolor="black")

    p2Y = []
    for n in np:
        my_lst = get_portfolio2(n)
        tts = []
        for l in my_lst:
            if l in cvc5_dict.keys():
                tts.append(cvc5_dict[l])
            else:
                tts.append(osmt_friendly[l])
        tt = get_total_time(tts)
        p2Y.append(sum(tt))

    plt.plot(np, p2Y, label="portfolio 2",
             alpha=.5, color="black", marker="o", linestyle="--", markerfacecolor="white", markeredgecolor="black")


    plt.xlabel("Number of Partitions")
    plt.ylabel("PAR-2 Score")
    plt.legend()
    plt.xscale("function", functions=(forward, inverse))
    plt.xticks(np)
    plt.savefig("figures/figure_4a.png", bbox_inches='tight')
    plt.clf()

    osmt_scatterY = [sum(osmt_dict[("scatter", n)]) for n in np]
    plt.plot(np, osmt_scatterY, label="OpenSMT2 scatter", linestyle="--", marker="X", color="black", alpha=.5)

    orpX = []
    orpY = []
    for n in np:
        my_lst = get_scatter_portfolio_osmt(n, osmt_dict)
        tt = get_total_time([osmt_dict[l] for l in my_lst])
        orpY.append(sum(tt))
        orpX.append(get_req_num_cores_osmt(my_lst))

    plt.plot(orpX, orpY, label="OpenSMT2 scatter portfolio (graduated)",
            alpha=.5, color="black", marker="*", markersize=10)
        
    plt.xlabel("Number of Partitions")
    plt.ylabel("PAR-2 Score")
    plt.legend()
    plt.xscale("function", functions=(forward, inverse))
    plt.xticks(np)
    plt.savefig("figures/figure_4b.png", bbox_inches='tight')
    plt.clf()

    osmt_scatterY = [sum(osmt_dict[("scatter", n)]) for n in np]
    plt.plot(np, osmt_scatterY, label="OpenSMT2 scatter", linestyle="--", marker="X", color="black", alpha=.5)

    orpX = []
    orpY = []
    for n in np:
        my_lst = get_scatter_portfolio_osmt(n, osmt_dict)
        tt = get_total_time([osmt_dict[l] for l in my_lst])
        orpY.append(sum(tt))
        orpX.append(get_req_num_cores_osmt(my_lst))

    plt.plot(orpX, orpY, label="OpenSMT2 scatter portfolio (graduated)",
             alpha=.5, color="black", marker="*", markersize=10)

    plt.plot(np, p2Y, label="portfolio 2",
             alpha=.5, color="black", marker="o", linestyle="--", markerfacecolor="white", markeredgecolor="black")

    PPX = [get_req_num_cores(get_graduated_portfolio_v1(n, cvc5_dict, osmt_friendly)) for n in np]
    PPY = []
    for n in np:
        lst = get_graduated_portfolio_v1(n, cvc5_dict, osmt_friendly)
        tt = []
        for l in lst:
            if l in cvc5_dict.keys():
                tt.append(cvc5_dict[l])
            else:
                tt.append(osmt_friendly[l])
        PPY.append(sum(get_total_time(tt)))


    plt.plot(PPX, PPY, label="portfolio 3 (graduated)", color="black", 
              marker="P", markeredgecolor="black", 
                markerfacecolor="white", alpha=.5, markersize=8)

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
        PPY.append(sum(get_total_time(tt)))

    plt.plot(PPX, PPY, label="recommended portfolio (graduated)", color="black", 
              marker="^", markeredgecolor="black", 
                markerfacecolor="white", alpha=.5, markersize=7)

    plt.xlabel("Number of Partitions")
    plt.ylabel("PAR-2 Score")
    plt.legend()
    plt.xscale("function", functions=(forward, inverse))
    plt.xticks([2,4,8,16,32,64,128])
    plt.savefig("figures/figure_4c.png", bbox_inches='tight')
    plt.clf()
    