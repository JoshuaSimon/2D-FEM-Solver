// FEM_GNUPlot.h
#include "pch.h"
#include "FEM_GNUPlot.h"

GNUPlot::GNUPlot()
{
	_pipe = nullptr;
}

GNUPlot::~GNUPlot()
{
	close();
}

void GNUPlot::open()
{
	close();
	_pipe = _popen("gnuplot", "w");
}

bool GNUPlot::isOpened() const
{
	return _pipe != nullptr;
}

void GNUPlot::flush()
{
	if (isOpened())
	{
		fflush(_pipe);
	}
}

void GNUPlot::close()
{
	if (isOpened())
	{
		_pclose(_pipe);
		_pipe = nullptr;
	}
}

void GNUPlot::write(const char *line)
{
	if (isOpened() && line != nullptr && line[0] != '\0')
	{
		fprintf(_pipe, "%s\n", line);
	}
}

void GNUPlot::write(const std::string &line)
{
	if (!line.empty())
	{
		write(line.c_str());
	}
}

void GNUPlot::execute(const std::vector<std::string> &script)
{
	if (isOpened())
	{
		for (size_t i = 0; i < script.size(); i++)
		{
			write(script[i]);
			flush();
		}
	}
}

//Function for plotting Mesh and Mesh + U
void plot(std::string filename)
{
	// Assemble plotting command
	std::string p_cmd;
	p_cmd.append("plot '");
	p_cmd.append(filename);
	p_cmd.append("' title \"Mesh\" with lines linetype 1 linecolor rgb \"#0000FF\" linewidth 3 , ");
	p_cmd.append(" '");
	p_cmd.append(filename);
	p_cmd.append("' using ($1 + 10 * $3) : ($2 + 10 * $4) title \"Mesh+U\" with lines linetype 1 \linecolor rgb \"#FF0000\" linewidth 3");

	// Commands send to gnuplot
	std::vector<std::string> script;
	script.push_back("set terminal wxt");
	script.push_back("reset");
	script.push_back("set xrange [-15:15]");
	script.push_back("set yrange [-15:15]");
	script.push_back("unset xtics");
	script.push_back("unset ytics");
	script.push_back(p_cmd);

	// Open gnuplot and plot graph
	GNUPlot plotter;
	plotter.open();
	plotter.execute(script);

	// Prevent graph to close
	system("Pause");

	// Close graph
	plotter.write("exit");
	plotter.flush();
	plotter.close();

	return;
}

