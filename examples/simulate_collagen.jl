using Plots
using SecondHarmonicGeneration
# using TestImages
# using Random

function main()

    # Random.seed!(1234)
    # Chi2 = rand(ComplexF64, 3,3,3)
    # Chi2 = rand(3,3,3)

    # chi2 = get_chi2_tensor(1)  # 3 for z-axis as fiber axis
    chi2 = get_chi2_stoller(1)  # 3 for z-axis as fiber axis
    αs = -pi:0.01:pi
    βs = 0 

    # surprisingly the chirality of the collagen seems to produce this non-symmetric?
    plot(αs.*180/pi, get_intensity.(Ref(chi2),αs, βs, false, true).*180/pi, xlabel="Angle / deg", yrange=(0,0.4),ylabel="Intensity", label="Collagen", title="Colinear Polarizer/Analyser")
    p2 = twinx()
    plot!(p2, αs.*180/pi, get_angle.(Ref(chi2),αs, βs, false).*180/pi, yrange=(0,180), ylabel="Phase Angle / deg", color=:red,  label=nothing)

    savefig("collagen_colinear.png")
    savefig("collagen_colinear.svg")

    plot(αs.*180/pi, get_intensity.(Ref(chi2),αs, βs, false, 1.0, false).*180/pi, xlabel="Angle / deg", yrange=(0,0.4), ylabel="Intensity", label="Collagen", title="Excitation Polarizer only")
    # p2 = twinx()
    # plot!(p2, αs.*180/pi, get_angle.(Ref(chi2),αs, βs, false).*180/pi, ylabel="Phase Angle / deg", color=:red,  label=nothing)

    savefig("collagen_expol.png")
    savefig("collagen_expol.svg")

    # circular polarizer/analyser
    cp = true
    axis_ratio = 1.0
    plot(αs.*180/pi, get_intensity.(Ref(chi2),αs, βs, cp, axis_ratio).*180/pi, xlabel="Angle / deg", yrange=(0,0.4), ylabel="Intensity", label="Collagen", title="Circular Polarizer/Analyser", legend=:topright)
    p2 = twinx()
    plot!(p2, αs.*180/pi, get_angle.(Ref(chi2),αs, βs, cp, axis_ratio).*180/pi, yrange=(-190,190), ylabel="Phase Angle / deg", color=:red, label=nothing)

    savefig("collagen_circular.png")
    savefig("collagen_circular.svg")

    cp = true
    axis_ratio = 0.4
    plot(αs.*180/pi, get_intensity.(Ref(chi2),αs, βs, cp, axis_ratio).*180/pi, xlabel="Angle / deg", yrange=(0,0.4), ylabel="Intensity", label="Collagen", title="Elliptic Polarizer/Analyser", legend=:topright)
    p2 = twinx()
    plot!(p2, αs.*180/pi, get_angle.(Ref(chi2),αs, βs, cp, axis_ratio).*180/pi, yrange=(-190,190), ylabel="Phase Angle / deg", color=:red, label=nothing)
    savefig("collagen_elliptic.png")
    savefig("collagen_elliptic.svg")

    plot(αs.*180/pi, hu_paper.(αs), xlabel="Angle / deg", ylabel="Intensity", label="Collagen", title="Hu paper")

    plot(αs.*180/pi, stoller_paper.(αs, pi/2), xlabel="Angle / deg", ylabel="Intensity", label="using Eqn. 8, γ", title="Stoller et al.")
    plot!(αs.*180/pi, 3.55*get_intensity.(Ref(chi2),αs, βs, false, false).*180/pi, xlabel="Angle / deg", ylabel="Intensity", label="using Eqn. 2, Χ²")

end
