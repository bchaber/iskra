using Circuit
v(t) = sin(2π*200*t)
cir = rlc(@netlist begin
	V1, VCC, GND, t -> 1.0
	L1, NOD, VCC, 500e-3
	C1, NOD, VCC, 3.5e-6
	R1, GND, NOD, 260
end)

data = zeros(500, 3)
Δt = 0.0001
for i=1:500
	data[i,:] = [cir.t, cir.i, cir.q]
	advance_circuit!(cir, 0.0, Δt)
end

using GR
plot(data[:,1], data[:,2])
print("press any key..."), readline()
