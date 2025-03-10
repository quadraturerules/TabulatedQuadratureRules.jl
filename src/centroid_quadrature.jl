function centroid_quadrature(
    domain::Domain,
    order::Integer,
)
    if domain == Domain_Interval && order == 1
        return [[0.5, 0.5]], [1.0]
    end
    if domain == Domain_Quadrilateral && order == 1
        return [[0.25, 0.25, 0.25, 0.25]], [1.0]
    end
    if domain == Domain_Triangle && order == 1
        return [[0.3333333333333333, 0.3333333333333333, 0.3333333333333333]], [1.0]
    end
    if domain == Domain_Hexahedron && order == 1
        return [[0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125]], [1.0]
    end
    if domain == Domain_Tetrahedron && order == 1
        return [[0.25, 0.25, 0.25, 0.25]], [1.0]
    end
    if domain == Domain_SquareBasedPyramid && order == 1
        return [[0.2, 0.2, 0.2, 0.2, 0.2]], [1.0]
    end
    if domain == Domain_TriangularPrism && order == 1
        return [[0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666]], [1.0]
    end
    throw("Invalid domain or order: ($domain, $order)")
end

export centroid_quadrature
