using DelimitedFiles

function lxcatread(filename, header)
	f = open(filename, "r")
	readuntil(f, header)
	h = readuntil(f, "-----------------------------")
	s = readuntil(f, "-----------------------------")
	readdlm(IOBuffer(s))
end
