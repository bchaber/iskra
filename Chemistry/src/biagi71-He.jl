import Chemistry: lxcatread, CrossSection
σ₁ = lxcatread("Biagi-7.1.txt", "PROCESS: E + He -> E + He, Elastic") |> CrossSection
σ₂ = lxcatread("Biagi-7.1.txt", "PROCESS: E + He -> E + He(Triplet)(19.82eV)") |> CrossSection
σ₃ = lxcatread("Biagi-7.1.txt", "PROCESS: E + He -> E + He(Singlet)(20.61eV)") |> CrossSection
σ₄ = lxcatread("Biagi-7.1.txt", "PROCESS: E + He -> E + E + He+, Ionization")  |> CrossSection
