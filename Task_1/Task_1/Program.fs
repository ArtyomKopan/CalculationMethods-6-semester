module Task1

open FsAlg.Generic
open FsUnit
open NUnit.Framework

let cond_v (A: Matrix<double>) (N: int) =
    let squares_sum n =
        List.fold (fun acc m -> acc + A.[n, m] ** 2) 0.0 [0..N-1]
        |> sqrt
    let product = List.fold (fun acc n -> acc * (squares_sum n)) 1.0 [0..N-1]
    product / (Matrix.det A)

let cond_a (A: Matrix<double>) (N: int) =
    let C = Matrix.inverse A
    List.map 
        (fun i -> Vector.norm (Matrix.toVector A[i,*]) * Vector.norm (Matrix.toVector C[*,i])) 
        [0..N-1]
    |> List.max 

let matrixNorm (a: Matrix<double>) (N: int) =
    let sumOfAbs j =
        List.fold (fun acc i -> acc + abs a.[i, j]) 0.0 [0..N-1]
    List.init N sumOfAbs |> List.max

let cond_s (A: Matrix<double>) (N: int) =
    (matrixNorm A N) * (matrixNorm (Matrix.inverse A) N)


let variate (A: Matrix<double>) varianceValue = A + varianceValue

let round3 (x: double) = System.Math.Round(x, 3)

let roundVector (v: Vector<double>) =
    Array.map round3 (Vector.toArray v) |> vector

[<Test>]
let test1() =
    let b = vector [200.0; -600.0]
    let A = matrix [[-400.6; 199.8]
                    [1198.8; -600.4]]
    
    printfn $"cond_v = %f{cond_v A 2}"
    printfn $"cond_a = %f{cond_a A 2}"
    printfn $"cond_s = %f{cond_s A 2}"
    
    let x = Matrix.solve A b |> roundVector
    let A' = variate A 1e-3
    let x' = Matrix.solve A' b |> roundVector
    printfn $"x = %A{x}"
    printfn $"x' = %A{x'}"
    printfn $"|x - x'| = %f{Vector.norm (x - x')}"
    
    round (cond_s A 2) |> should equal 2878


[<Test>]
let test2() =
    // матрица Гильберта
    let A = matrix [[1.0; 0.5; 1.0/3.0]
                    [0.5; 1.0/3.0; 0.25]
                    [1.0/3.0; 0.25; 0.2]]
    let b = vector [1.0; 1.0; 1.0]
    let A' = variate A 1e-3
    let x = Matrix.solve A b |> roundVector
    let x' = Matrix.solve A' b |> roundVector
    
    printfn $"cond_v = %f{cond_v A 3}"
    printfn $"cond_a = %f{cond_a A 3}"
    printfn $"cond_s = %f{cond_s A 3}"
    printfn $"x = %A{x}"
    printfn $"x' = %A{x'}"
    printfn $"|x - x'| = %f{Vector.norm (x - x')}"
    
    round (cond_s A 3) |> should equal 748
    
[<Test>]
let test3() =
    // трёхдиагональная матрица с диагональным преобладанием
    let A = matrix [[6.; 3.; 0.; 0.;0.]
                    [1.; 5.; 2.; 0.; 0.;]
                    [0.; 4.; 10.; 3.; 0.]
                    [0.; 0.; 2.; 5.; 1.]
                    [0.; 0.; 0.; 3.; 4.]]
    let b = vector [1.; 1.; 1.; 1.; 1.]
    let A' = variate A 1e-3
    let x = Matrix.solve A b |> roundVector
    let x' = Matrix.solve A' b |> roundVector
    
    printfn $"cond_v = %f{cond_v A 5}"
    printfn $"cond_a = %f{cond_a A 5}"
    printfn $"cond_s = %f{cond_s A 5}"
    printfn $"x = %A{x}"
    printfn $"x' = %A{x'}"
    printfn $"|x - x'| = %f{Vector.norm (x - x')}"
    
    round (cond_s A 5) |> should equal 9