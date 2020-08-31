import Chemistry: lxcatread, CrossSection
σₑ₁ = lxcatread("Phelps.txt", "PROCESS: He+ + He -> , Isotropic") |> Chemistry.CrossSection
σₑ₂ = lxcatread("Phelps.txt", "PROCESS: He+ + He -> , Backscat")  |> Chemistry.CrossSection	
