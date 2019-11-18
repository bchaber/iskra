using Circuit
v(t) = sin(2π*200*t)
cir = rlc(@netlist begin
	V1, VCC, GND, 3.0
	R1, GND, NOD, 10 + 40
	L1, NOD, VCC, 1e-6
end)