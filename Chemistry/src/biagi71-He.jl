import Chemistry: lxcatread, CrossSection
σ₁ = lxcatread("datasets/Biagi-7.1.txt", "PROCESS: E + He -> E + He, Elastic") |> CrossSection
σ₂ = lxcatread("datasets/Biagi-7.1.txt", "PROCESS: E + He -> E + He(Triplet)(19.82eV)") |> CrossSection
σ₃ = lxcatread("datasets/Biagi-7.1.txt", "PROCESS: E + He -> E + He(Singlet)(20.61eV)") |> CrossSection
σ₄ = lxcatread("datasets/Biagi-7.1.txt", "PROCESS: E + He -> E + E + He+, Ionization")  |> CrossSection
