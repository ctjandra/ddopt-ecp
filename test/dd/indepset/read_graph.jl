using LightGraphs

"""
Read a graph in DIMACS format. Assumes file is correctly formatted.
"""
function read_graph(filename::AbstractString)::Graph
	f = open(filename)
	g = nothing
	for line in eachline(f)
		pieces = split(line)			#split strings into pieces separated by space (default delimiter)
		if pieces[1] == "p"
			nvertices = parse(Int, pieces[3])		#converts string to int values
			g = Graph(nvertices)
		elseif pieces[1] == "e"
			if g == nothing
				error("clq file does not define number of vertices in header")
			end
			add_edge!(g, parse(Int, pieces[2]), parse(Int, pieces[3]))
		end
	end
	return g
end
