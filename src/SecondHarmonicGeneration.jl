module SecondHarmonicGeneration


export get_chi2_tensor, get_chi2_stoller, get_EP, get_intensity, get_angle, hu_paper, stoller_paper

using Tullio
using LinearAlgebra

"""
Gemini:
For a material with χ²_{∞v} symmetry, if we align the z-axis with the fiber axis (the axis of the triple helix/fibril),
the non-zero independent elements of the χ^{(2)} tensor are:
 * χ²_{zzz}
 * χ²_{zxx} (where x and y are arbitrary directions perpendicular to z)
 * χ²_{zyy} (which is equal to χ²_{zxx} under C_{∞v} symmetry)
 * χ²_{xxz} (which is equal to χ²_{xzz}, χ²_{yyz}, χ²_{yzz} due to Kleinman symmetry and C_{∞v} symmetry)
Due to Kleinman symmetry (which applies when there is no absorption at any of the involved frequencies,
a good approximation at 800nm/400nm for collagen), the tensor is further simplified.
Kleinman symmetry states that the permutation of all indices (including frequency arguments) does not change
the tensor elements. This implies:
χ²_{ijk} = χ²_{jik} = χ²_{ikj} etc.
Applying Kleinman symmetry to the C_{∞v} non-zero elements, we get the following relationships:
 * χ²_{zxx} = χ²_{xxz} (and similarly for y)
So, for collagen with C_{∞v} symmetry and assuming Kleinman symmetry, there are generally two independent non-zero coefficients:
 * χ²_{zzz}: Corresponds to the fundamental electric field polarized along the collagen fiber axis (z) generating SHG polarized
 along the fiber axis.
 * χ²_{zxx} (or χ²_{zyy}): Corresponds to the fundamental electric field polarized perpendicular to the fiber axis (x or y) generating SHG polarized along the fiber axis, OR the fundamental electric field having a component along z and a component along x (or y) generating SHG along x (or y).
Often, researchers also talk about the ratio ρ = χ²_{zzz} / χ²_{zxx} as a key parameter for characterizing collagen.
This ratio typically ranges from ~0.5 to ~2.0 depending on the collagen type and structural packing.
"""
function get_chi2_tensor(orientation_dim=1)
    χ² = zeros(ComplexF64, 3, 3, 3)

    # type 1 collagen as mostly found in dura mater
    axis_value = 0.07 # pV/m
    off_axis_value = axis_value / 1.4 # pV/m   for zxx
    other_off_axis_value = off_axis_value # pV/m  ??
    # The ratio should be between 1.5 and 2.5
    # first construct the tensor with z (orientation_dim=3) as the fiber axis

    z = 3; x = 1; y = 2
    if orientation_dim == 1
        z = 1; x = 3; y = 2
    elseif orientation_dim == 2
        z = 2; x = 1; y = 3
    elseif orientation_dim != 3
        error("orientation_dim must be 1, 2, or 3")
    end 

    # assuming Kleinman symmetry and C_{∞v} symmetry:
    χ²[z, z, z] = axis_value 

    χ²[z, x, x] = off_axis_value 
    χ²[z, y, y] = off_axis_value 
    χ²[x, z, x] = off_axis_value 
    χ²[y, z, y] = off_axis_value 
    χ²[x, x, z] = off_axis_value 
    χ²[y, y, z] = off_axis_value 

    # χ²[x, z, z] = off_axis_value 
    # χ²[y, z, z] = off_axis_value
    # χ²[z, x, z] = off_axis_value 
    # χ²[z, y, z] = off_axis_value 
    # χ²[z, z, x] = off_axis_value 
    # χ²[z, z, y] = off_axis_value 

    # if the orientation is not z, we need to permute the tensor
    return χ²
end

"""
    get_chi2_stoller(orientation_dim=1, gamma=-0.7, a=0.07, b = gamma*a, c=2*b)

Get the χ² tensor according equation (2) of Stoller et al. (2002) for collagen type 1.

gamma = -0.7 = b/a
b = gamma * a
b = c/2
"""
function get_chi2_stoller(orientation_dim=1, gamma=-0.7, a=0.07, b = gamma*a, c=2*b)
    # Kleinman symmetry warrants that b = c/ 2 according to Stoller et al.
    if (orientation_dim == 1)
        s = [1,0,0]
    elseif (orientation_dim == 2)
        s = [0,1,0]
    elseif (orientation_dim == 3)
        s = [0,0,1]
    else
        error("orientation_dim must be 1, 2, or 3")
    end

    χ² = zeros(ComplexF64, 3, 3, 3)
    @tullio χ²[i,j,k] = a*s[i]*s[j]*s[k]  + b * s[i] *(j==k) + c/2*(s[j]*(i==k) + s[k]*(i==j))  # pV/m

    return χ²
end

function get_EP(chi2, α=0, β=0, circular_polarization=false, axis_ratio=1.0)
    E = [cos(α)*cos(β), sin(α)*cos(β), sin(β)]
    if (circular_polarization)
        # circular polarization
        E = (1/sqrt(1+abs2(axis_ratio))) .*(E .+ axis_ratio.*1im.*[sin(α)*cos(β), cos(α)*cos(β), sin(β)])  # convert to circular polarization
    end
    P = zeros(ComplexF64, 3)
    @tullio P[k] = chi2[k,i,j] * E[i] * E[j]
    return E,P
end

function get_intensity(chi2, α=0, β=0, circular_polarization=false,  axis_ratio=1.0, analyzer = true)
    E,P = get_EP(chi2, α, β, circular_polarization, axis_ratio)

    # E_det = [cos(α)*cos(β), sin(α)*cos(β), sin(β)]
    E_det = E
    if (analyzer)
        return abs2(sum(P .* E_det))
    else
        return sum(abs2.(P))
    end
end

function get_angle(chi2, α=0, β=0, circular_polarization=false,  axis_ratio=1.0)
    E,P = get_EP(chi2, α, β, circular_polarization, axis_ratio)
    # E_det = [cos(α)*cos(β), sin(α)*cos(β), sin(β)]
    E_det = E
    measured = sum(P .* E_det)
    return angle(measured)
end

function hu_paper(phi)
    d33  = 0.07
    d31 = d33 / 1.4
    d15 = d31 * 0.53

    Px = d15 * sin(2*phi)
    Pz = d31 * abs2(sin(phi)) + d33 * abs2(cos(phi))

    I = abs2(Px) + abs2(Pz)
    return I
end

"""
    stoller_paper(alpha, phi=0, gamma=-0.7)

this uses a different equation (8) from the Stoller paper, but it should be equivalent.

"""
function stoller_paper(alpha, phi=0, gamma=-0.7)
    I = 1/8 * (3 + 20*gamma + 40*gamma^2) -
        1/2 * (1 + 6*gamma +8*gamma^2)*cos(2*alpha + 2*phi) + 
        1/8*(1 +4*gamma)*cos(4*alpha + 4*phi)

    return I
end

end # module SecondHarmonicGeneration
