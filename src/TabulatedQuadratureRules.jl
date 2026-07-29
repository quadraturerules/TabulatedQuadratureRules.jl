module TabulatedQuadratureRules

include("domains.jl")
include("clenshaw_curtis.jl")
include("gauss_legendre.jl")
include("gauss_lobatto_legendre.jl")
include("hammer_marlowe_stroud.jl")
include("sauter_schwab.jl")
include("xiao_gimbutas.jl")
include("centroid_quadrature.jl")
include("closed_newton_cotes.jl")
include("open_newton_cotes.jl")
include("vertex_quadrature.jl")

@enum QuadratureRule begin
    QR_ClenshawCurtis = 10
    QR_GaussLegendre = 1
    QR_GaussLobattoLegendre = 3
    QR_HammerMarloweStroud = 6
    QR_SauterSchwab = 7
    QR_XiaoGimbutas = 2
    QR_CentroidQuadrature = 5
    QR_ClosedNewtonCotes = 8
    QR_OpenNewtonCotes = 9
    QR_VertexQuadrature = 4
end

"Get a quadrature rule for a single integral."
function single_integral_quadrature(
    rtype::QuadratureRule,
    domain::Domain,
    order::Integer,
)
    if rtype == QR_ClenshawCurtis
        return clenshaw_curtis(domain, order)
    end
    if rtype == QR_GaussLegendre
        return gauss_legendre(domain, order)
    end
    if rtype == QR_GaussLobattoLegendre
        return gauss_lobatto_legendre(domain, order)
    end
    if rtype == QR_HammerMarloweStroud
        return hammer_marlowe_stroud(domain, order)
    end
    if rtype == QR_XiaoGimbutas
        return xiao_gimbutas(domain, order)
    end
    if rtype == QR_CentroidQuadrature
        return centroid_quadrature(domain, order)
    end
    if rtype == QR_ClosedNewtonCotes
        return closed_newton_cotes(domain, order)
    end
    if rtype == QR_OpenNewtonCotes
        return open_newton_cotes(domain, order)
    end
    if rtype == QR_VertexQuadrature
        return vertex_quadrature(domain, order)
    end
    throw("Unsupported rule for single integral: $rtype")
end

"Get a quadrature rule for a double integral."
function double_integral_quadrature(
    rtype::QuadratureRule,
    domain::Domain,
    order::Integer,
)
    if rtype == QR_SauterSchwab
        return sauter_schwab(domain, order)
    end
    throw("Unsupported rule for single integral: $rtype")
end

export QuadratureRule
export Domain
export single_integral_quadrature

end
