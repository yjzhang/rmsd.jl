
# function for loading coordinate data
read_structure(filename) = float64(readdlm(filename)[2:end,:])

function write_structure(coords::Array{Float64, 2}, out_file::String)
    # coords: each row is a coord
    # outputs a space-separated text file
    n, dim = size(coords)
    out = open(out_file, "w")
    write(out, @sprintf("%d\n", n))
    for i in 1:n
        c = coords[i,:]
        write(out, @sprintf("%8.6f ", c[1]))
        write(out, @sprintf("%8.6f ", c[2]))
        write(out, @sprintf("%8.6f ", c[3]))
        write(out, "\n")
    end
    close(out)
end

function rmsd_align(v1::Array{Float64, 2}, v2::Array{Float64, 2})
    # returns copy of v1 rotated to align to v2
    (_, R, _) = rmsd_rotation(v1, v2)
    rotated_v1 = (R*v1)'
    return rotated_v1
end


function rmsd_rotation(v1::Array{Float64, 2}, v2::Array{Float64, 2})
    # v1 and v2 are N x 3 arrays
    # basically, this is the Kabsch algorithm (see Wiki)
    # returns the RMSD, rotation matrix, and translation
    c1 = mean(v1,1)
    c2 = mean(v2,1)
    v1_centered = v1.-c1
    v2_centered = v2.-c2
    A = (v1_centered)'*(v2_centered)
    (U, S, V) = svd(A)
#    println(U)
#    println(S)
#    println(V)
    R = V*U'
#    if det(R) < 0
#        I = eye(3)
#        I[3,3] = -1.0
#        R = V*I*U'
#    end
    v1_rotated = R*v1'
    translation = -R*c1' + c2'
    v1_final = v1_rotated .+ translation
    return rmsd(v1_final', v2), R, translation
end

function rmsd_horn(v1::Array{Float64,2}, v2::Array{Float64,2})
    # Calculates optimal rotation, scaling, translation based on Horn (1988)
    c1 = mean(v1,1)
    c2 = mean(v2,1)
    v1_centered = v1.-c1
    v2_centered = v2.-c2

    M = (v1_centered)'*(v2_centered)
    M = M'
    MTM = M'*M
    eigval, eigvec = eig(MTM)
    S_inv = 1/sqrt(eigval[1])*eigvec[:,1]*eigvec[:,1]' + 
            1/sqrt(eigval[2])*eigvec[:,2]*eigvec[:,2]' +
            1/sqrt(eigval[3])*eigvec[:,3]*eigvec[:,3]'
    R = M*S_inv
    s = sqrt(sum(sum(v1_centered, 2).^2)/sum(sum(v2_centered,2).^2))
    #s = 1
    t = -s*R*c1' + c2'
    v1_final = s*R*v1' .+ t
    return rmsd(v1_final', v2), R, s, t
end

function distance_matrix(points::Array{Float64,2})
    # Returns a distance matrix given an array of points N x 3.
    dist = zeros(size(points,1),size(points,1))
    for i = 1:size(points,1)
        for j = (i+1):size(points,1)
            dist[i,j] = sqrt(sum((points[i,:]-points[j,:]).^2))
            dist[j,i] = dist[i,j]
        end
    end
    return dist
end


function rmsd(v1::Array{Float64, 2}, v2::Array{Float64, 2})
    # returns the RMSD between two arrays, where each row
    # is a point.
    return sqrt(sum(sum(v1 - v2, 2).^2))/size(v1,1)
end
