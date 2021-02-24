import Chemistry: lxcatread, CrossSection
σ₁ = lxcatread("datasets/Biagi-7.1.txt", "PROCESS: E + Ar -> E + Ar, Elastic") |> CrossSection
σ₂ = lxcatread("datasets/Biagi-7.1.txt", "PROCESS: E + Ar -> E + Ar*(11.55eV), Excitation") |> CrossSection
σ₃ = lxcatread("datasets/Biagi-7.1.txt", "PROCESS: E + Ar -> E + Ar*(13.00eV), Excitation") |> CrossSection
σ₄ = lxcatread("datasets/Biagi-7.1.txt", "PROCESS: E + Ar -> E + E + Ar+, Ionization")  |> CrossSection
