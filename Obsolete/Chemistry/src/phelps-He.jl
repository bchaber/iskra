import Chemistry: lxcatread, CrossSection
σₑ₁ = lxcatread("datasets/Phelps.txt", "PROCESS: He+ + He -> , Isotropic") |> Chemistry.CrossSection
σₑ₂ = lxcatread("datasets/Phelps.txt", "PROCESS: He+ + He -> , Backscat")  |> Chemistry.CrossSection	
