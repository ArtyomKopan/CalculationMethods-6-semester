open FsAlg.Generic

let cond_v (A: Matrix<double>) (N: int) =
    let squares_sum n =
        sqrt (List.fold (fun acc m -> acc + A.[n, m] ** 2) 0.0 [0..N-1])
    let product = List.fold (fun acc n -> acc * (squares_sum n)) 1.0 [0..N-1]
    product / (Matrix.det A)

let cond_a (A: Matrix<double>) (N: int) =
    let C = Matrix.inverse A
    List.max (List.map 
        (fun i -> Vector.norm (Matrix.toVector A[i,*]) * Vector.norm (Matrix.toVector C[*,i])) 
        [0..N-1])

let matrixNorm (a: Matrix<double>) (N: int) =
    let sumOfAbs j =
        List.fold (fun acc i -> acc + abs a.[i, j]) 0.0 [0..N-1]
    List.fold (fun acc j -> acc * sumOfAbs j) 1.0 [0..N-1]

let cond_s (A: Matrix<double>) (N: int) =
    (matrixNorm A N) * (matrixNorm (Matrix.inverse A) N)


let A = matrix [[1.0; -1.0]; [2.0; 1.0]]
let b = vector [-5.0; -7.0]

let x = Matrix.solve A b

printfn "%A" x