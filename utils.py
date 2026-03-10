import numpy as np
from scipy.stats import norm, hypergeom, fisher_exact
from statsmodels.stats.proportion import proportions_ztest
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import seaborn as sns
prop = font_manager.FontProperties(fname = "/home/tx82/Arial.ttf")

# --------------------------------------------------------------
# For stats
# --------------------------------------------------------------
def calc_OR(k, N, n, M, alternative="greater"):
	"""
	Calculates odds ratio, confidence interval, and p value.
	Perform Fisher's exact test.

	Parameters
	----------
	k: integer
		The number of observed successes (e.g. intersection between candidate set and benchmark set).
	N: integer
		The number of observed draws (e.g. candidate gene set).
	n: integer
		The number of all successes in the population (e.g. benchmark gene set).
	M: integer
		The population size (e.g. all genes).

	Returns
	-------
	OR: float
		odds ratio, np.nan if invalid
	CI: (float, float)
		95% confidence interval for the odds ratio, np.nan if invalid
	pval: float
		p-value for the null hypothesis that the log odds ratio is 0.
	"""
	OR, pval = fisher_exact([[k, n-k], [N-k, M-(n+N)+k]], alternative = alternative)
	try:
		lower_ci = np.exp(np.log(OR) - 1.96 * np.sqrt(1 / k + 1 / (n-k) + 1 / (N-k) + 1 / (M-(n+N)+k)))
	except:
		lower_ci = np.nan
	try:
		upper_ci = np.exp(np.log(OR) + 1.96 * np.sqrt(1 / k + 1 / (n-k) + 1 / (N-k) + 1 / (M-(n+N)+k)))
	except:
		upper_ci = np.nan
	CI = (lower_ci, upper_ci)
	return [OR, CI, pval]

def proportion2test(x1, n1, x2, n2, alternative="two-sided"):
	'''
	Test whether the percentages in two samples are the same; alternative can be "larger"/"smaller"/"two-sided"
	'''
	import numpy as np
	if n1 == 0 or n2 == 0: 
		return np.nan
	else:
		_, p = proportions_ztest([x1, x2], [n1, n2], alternative=alternative)
		return p

# --------------------------------------------------------------
# For plot
# --------------------------------------------------------------
def format_pval(pval, pval_star = False, pval_cutoff = 0.05, pval_sci = True):
    '''
    Format p value for plotting
    '''
    if pval_star == False:
        if pval >= pval_cutoff or pval is None:
            label = "n.s."
        else:
            if pval_sci == True:
                label = '{:0.1e}'.format(pval)
            else:
                label = str(round(pval, 2))
    else:
        if pval >= pval_cutoff:
            label = "n.s."
        elif pval is None or np.isnan(pval):
            label = None
        else:
            if pval == 0:
                label = "*" * 3
            elif (int('{:0.1e}'.format(pval).split("-")[1]) - 1) <= 3:
                label = "*" * (int('{:0.1e}'.format(pval).split("-")[1]) - 1)
            else:
                label = "*" * 3

    return label

def single_bar_graph(df, x_col, y_col, pval_col=None, pval_fontsize=10, pval_parameters={"pval_cutoff":0.05, "pval_sci":False, "pval_star":True},\
                     y_err_col=None, y_lerr_col=None, y_uerr_col=None, output=None, \
                      bar_label_col=None, bar_label_parameters={"label_type":"edge", "fontsize":8, "fontweight":"semibold", "padding":1}, \
                      figsize=None, fontsize=11, color=sns.color_palette(palette="Set2")[2], palette=None, legend_parameters={"fontsize":11, "loc":"best", "bbox_to_anchor":(1, 0.5), "handlelength":1}, \
                        xlim=None, ylim=None, xticks_rotation=None, yticks_rotation=None, xtick_align_ha="center", xtick_align_va="top", xlabel=None, ylabel=None, title=None, ax=None):
	'''
	Plot a single bar graph with optional error bars, bar labels, and p-value labels
	'''
	prop.set_size(fontsize)
	# prop["size"] = fontsize

	fig = plt.figure(figsize=figsize)
	if ax is not None:
		sns.barplot(x=x_col, y=y_col, data=df, color=color, palette=palette, ax=ax)
	else: 
		ax = sns.barplot(x=x_col, y=y_col, data=df, color=color, palette=palette)
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
	if legend_parameters is not None:
		ax.legend(prop=font_manager.FontProperties(fname = "/home/tx82/Arial.ttf"), **legend_parameters)
	ax.set_xticklabels([i.get_text() for i in ax.get_xticklabels()], rotation=xticks_rotation, ha=xtick_align_ha, va=xtick_align_va, fontproperties = prop)
	ax.set_yticklabels([i.get_text() for i in ax.get_yticklabels()], rotation=yticks_rotation, fontproperties = prop)
	if xlabel is not None:
		ax.set_xlabel(xlabel, fontproperties = prop)
	if ylabel is not None:
		ax.set_ylabel(ylabel, fontproperties = prop)
	if title is not None:
		ax.set_title(title, fontproperties = prop)
	
	# plot bar labels
	if bar_label_col is not None:
		bar_label_list = list(df[bar_label_col])
		for i,container in enumerate(ax.containers):
			if isinstance(container, matplotlib.container.BarContainer):
				ax.bar_label(container, labels=[bar_label_list[i]], **bar_label_parameters)

	# plot error bars
	if y_err_col is not None:
		y_err_list = list(df[y_err_col])
		for i,container in enumerate(ax.containers):
			if isinstance(container, matplotlib.container.BarContainer):
				x_coords = [p.get_x() + 0.5 * p.get_width() for p in container.patches]
				y_coords = [p.get_height() for p in container.patches]
				ax.errorbar(x=x_coords, y=y_coords, yerr=y_err_list[i], fmt="none", c="k")
	elif y_lerr_col is not None and y_uerr_col is not None:
		y_lerr_list = list(df[y_lerr_col])
		y_uerr_list = list(df[y_uerr_col])
		for i,container in enumerate(ax.containers):
			if isinstance(container, matplotlib.container.BarContainer):
				x_coords = [p.get_x() + 0.5 * p.get_width() for p in container.patches]
				y_coords = [p.get_height() for p in container.patches]
				ax.errorbar(x=x_coords, y=y_coords, yerr=np.array([y_lerr_list[i], y_uerr_list[i]]), fmt="none", c="k")

	# plot p-values
	if pval_col is not None:
		pval_labels = [format_pval(i, **pval_parameters) for i in list(df[pval_col])]
		for i,container in enumerate(ax.containers):
			if isinstance(container, matplotlib.container.BarContainer):
				y_coords = [p.get_height() for p in container.patches]
				if y_err_col is not None: 
					y_coords = list(np.array(y_coords) + np.array(list(df[y_err_col])))
				elif y_lerr_col is not None and y_uerr_col is not None:
					y_coords = list(np.array(y_coords) + np.array(list(df[y_uerr_col])))
				
				gap = (ax.get_ylim()[1] - ax.get_ylim()[0])*0.01
				for bar,label,y_coord in zip(container, [pval_labels[i]], y_coords):
						x_coord = bar.get_x() + 0.5 * bar.get_width()
						ax.text(x_coord, y_coord+gap, label, ha="center", fontsize=pval_fontsize, fontstyle="italic")
	
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	if output is not None:
		plt.savefig(output, bbox_inches='tight', dpi=500)
	
	return fig, ax

def grouped_bar_graph(df, x_col, y_col, hue_col, hue_order=None, pval_col=None, pval_fontsize=10, pval_parameters={"pval_cutoff":0.05, "pval_sci":False, "pval_star":True},\
                     y_err_col=None, y_lerr_col=None, y_uerr_col=None, output=None, \
                      bar_label_col=None, bar_label_parameters={"label_type":"edge", "fontsize":8, "fontweight":"semibold", "padding":1}, \
                      figsize=None, fontsize=11, palette_list=sns.color_palette(palette="Set2"), legend_parameters={"fontsize":11, "loc":"best", "bbox_to_anchor":(1, 0.5), "handlelength":1}, \
                        xlim=None, ylim=None, xticks_rotation=None, yticks_rotation=None, xtick_align_ha="center", xtick_align_va="top", xlabel=None, ylabel=None, title=None, ax=None):
      
	'''
    Generates a grouped bar plot with optional error bars, labels, and p-values
    '''
	
	prop.set_size(fontsize)
	# prop["size"] = fontsize

	fig = plt.figure(figsize=figsize)
	if ax is not None:
		sns.barplot(x=x_col, y=y_col, hue=hue_col, hue_order=hue_order, data=df, palette=palette_list, ax=ax)
	else: 
		ax = sns.barplot(x=x_col, y=y_col, hue=hue_col, hue_order=hue_order, data=df, palette=palette_list)
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
	prop_legend = font_manager.FontProperties(fname = "/home/tx82/Arial.ttf")
	if "fontsize" in legend_parameters:
		prop_legend.set_size(legend_parameters["fontsize"])
	ax.legend(prop=prop_legend, **legend_parameters)
	ax.set_xticklabels([i.get_text() for i in ax.get_xticklabels()], rotation=xticks_rotation, ha=xtick_align_ha, va=xtick_align_va, fontproperties = prop)
	ax.set_yticklabels([i.get_text() for i in ax.get_yticklabels()], rotation=yticks_rotation, fontproperties = prop)
	if xlabel is not None:
		ax.set_xlabel(xlabel, fontproperties = prop)
	if ylabel is not None:
		ax.set_ylabel(ylabel, fontproperties = prop)
	if title is not None:
		ax.set_title(title, fontproperties = prop)
	
	# plot bar labels
	if bar_label_col is not None:
		for container, group in zip(ax.containers, df[hue_col].unique()):
			if isinstance(container, matplotlib.container.BarContainer):
				df_container = df[df[hue_col] == group]
				ax.bar_label(container, labels=list(df_container[bar_label_col]), **bar_label_parameters)

	# plot error bars
	if y_err_col is not None:
		for container, group in zip(ax.containers, df[hue_col].unique()):
			if isinstance(container, matplotlib.container.BarContainer):
				df_container = df[df[hue_col] == group]
				x_coords = [p.get_x() + 0.5 * p.get_width() for p in container.patches]
				y_coords = [p.get_height() for p in container.patches]
				y_err_list = list(df_container[y_err_col])
				ax.errorbar(x=x_coords, y=y_coords, yerr=y_err_list, fmt="none", c="k")
	elif y_lerr_col is not None and y_uerr_col is not None:
		for container, group in zip(ax.containers, df[hue_col].unique()):
			if isinstance(container, matplotlib.container.BarContainer):
				df_container = df[df[hue_col] == group]
				x_coords = [p.get_x() + 0.5 * p.get_width() for p in container.patches]
				y_coords = [p.get_height() for p in container.patches]
				y_lerr_list = list(df_container[y_lerr_col])
				y_uerr_list = list(df_container[y_uerr_col])
				ax.errorbar(x=x_coords, y=y_coords, yerr=np.array([y_lerr_list, y_uerr_list]), fmt="none", c="k")

	# plot p-values
	if pval_col is not None:
		for container, group in zip(ax.containers, df[hue_col].unique()):
			if isinstance(container, matplotlib.container.BarContainer):
				df_container = df[df[hue_col] == group]
				y_coords = [p.get_height() for p in container.patches]
				if y_err_col is not None: 
					y_coords = list(np.array(y_coords) + np.array(list(df_container[y_err_col])))
				elif y_lerr_col is not None and y_uerr_col is not None:
					y_coords = list(np.array(y_coords) + np.array(list(df_container[y_uerr_col])))
				pval_labels = [format_pval(i, **pval_parameters) for i in list(df_container[pval_col])]

				for bar,label,y_coord in zip(container, pval_labels, y_coords):
						x_coord = bar.get_x() + 0.5 * bar.get_width()
						ax.text(x_coord, y_coord, label, ha="center", fontsize=pval_fontsize, fontstyle="italic")
	
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	if output is not None:
		plt.savefig(output, bbox_inches='tight', dpi=500)
	
	return fig, ax

def grouped_scatter_with_se(df, x_col, y_col, hue_col, hue_order, yerr_col=None, figsize=(8,3), fontsize=14,
    face_palette_list=None, edge_palette_list=None, xlabel="", ylabel="", title="", ylim=None, legend_parameters=None):
    '''
    Plot a grouped scatter graph with optional standard error bars
	'''
    # category order
    regions = df[x_col].drop_duplicates().tolist()
    n_regions = len(regions)

    x = np.arange(n_regions)
    n_hue = len(hue_order)
    total_width = 0.8
    box_width = total_width / n_hue

    offsets = {h: x - total_width/2 + (i + 0.5)*box_width for i, h in enumerate(hue_order)}

    if face_palette_list is None:
        face_palette_list = sns.color_palette("Pastel1", n_colors=n_hue)
    face_color_map = {h: face_palette_list[i] for i, h in enumerate(hue_order)}

    if edge_palette_list is None:
        edge_palette_list = sns.color_palette("Dark2", n_colors=n_hue)
    edge_color_map = {h: edge_palette_list[i] for i, h in enumerate(hue_order)}

    fig, ax = plt.subplots(figsize=figsize)

    # plot points + error bars
    for h in hue_order:
        sub = df[df[hue_col] == h].set_index(x_col).reindex(regions)
        y = sub[y_col].values.astype(float)
        yerr = sub[yerr_col].values.astype(float)

        ax.errorbar(
            offsets[h], y, yerr=yerr,
            fmt='o', markersize=7,
            color=face_color_map[h], markeredgecolor=edge_color_map[h], ecolor=edge_color_map[h],
            elinewidth=1.2, capsize=3, capthick=1.2,
            label=h, zorder=3
        )

    # x ticks
    ax.set_xticks(x)
    ax.set_xticklabels(regions, fontsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)

    # grid helps readability for points
    ax.yaxis.grid(True, linestyle='-', alpha=0.25)
    ax.set_axisbelow(True)

    if ylim is not None:
        ax.set_ylim(ylim)

    # legend
    if legend_parameters is None:
        legend_parameters = {"fontsize":12, "loc":"best", "bbox_to_anchor":(1,0.75), "handlelength":1}
    ax.legend(**legend_parameters)

    return fig, ax, offsets, regions