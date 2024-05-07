open FsAlg.Generic
open Plotly.NET
open System

type Circle(centerX: float, centerY: float, radius: float) =
    member c.centerX = centerX
    member c.centerY = centerY
    member c.radius = radius

let createHilbertMatrix n =
    Matrix.init n n (fun i j -> 1.0 / (float(i) + float(j) + 1.0))

let matrixFromList inputList =
    let numRows = List.length inputList
    let numCols = List.length (List.head inputList)
    let array2D = Array2D.init numRows numCols (fun i j -> inputList.[i].[j])
    Matrix.ofArray2D array2D

let GershgorinCircles (A: Matrix<float>) =
    let n = A.Rows
    let R = List.map (fun i -> 
                      (Matrix.row i A 
                      |> Matrix.toVector
                      |> Vector.toArray 
                      |> Array.map abs
                      |> Array.sum) - abs A[i, i]
                      ) [0..n - 1]
    let circles = List.map (fun i -> Circle(A[i, i], 0.0, R[i])) [0..n - 1]
    circles

let GershgorinCirclesCheck (A: Matrix<float>) (eigenValues: List<float>) =
    let n = A.Rows
    let circles = GershgorinCircles A
    [0..n - 1]
    |> List.forall (fun i -> abs (eigenValues[i] - circles[i].centerX) < circles[i].radius)

let chooseMaxAbsElement (A: Matrix<float>) =
    let n = A.Rows
    let mutable maxI = 0
    let mutable maxJ = 1
    for i in 0..n - 1 do
        for j in 0..n - 1 do
            if i <> j && abs A[i, j] > abs A[maxI, maxJ] then
                maxI <- i
                maxJ <- j
        done
    done
    maxI, maxJ

let chooseByMaxGershgorinCircle (A: Matrix<float>) =
    let n = A.Rows
    let circles = GershgorinCircles A
    let maxI = [0..n - 1]
               |> List.maxBy (fun i -> circles[i].radius)
    let mutable maxJ = 0
    for j in 1..n - 1 do
        if j <> maxI then
            if A[maxI, j] > A[maxI, maxJ] then
                maxJ <- j
    done 
    maxI, maxJ


let JacobiMethod (B: Matrix<float>) chooseStrategy epsilon maxIterations =
    
    let n = B.Rows

    let nondiagonalSquaresSum (A: Matrix<float>) =
        [0..n - 1]
        |> List.map (fun i -> 
                    (Matrix.row i A 
                    |> Matrix.toVector
                    |> Vector.toArray
                    |> Array.fold (fun acc elem -> acc + elem ** 2) 0.0
                    ) - A[i, i] ** 2)
        |> List.sum

    let rec step (A: Matrix<float>) numIteration =
        let (i, j) = chooseStrategy A
        let x = -2.0 * A[i, j]
        let y = A[i, i] - A[j, j]
        let c = if (abs y < 1e-6) then 1.0 / (sqrt 2.0) 
                  else 0.5 * (1.0 + (abs y) / (sqrt (x ** 2 + y ** 2))) |> sqrt
        let s = if (abs y < 1e-6) then 1.0 / (sqrt 2.0)
                  else (float (Math.Sign (x * y))) * (abs x) / (2.0 * c * sqrt (x ** 2 + y ** 2))
        let T = Matrix.init n n (fun i j -> if i = j then 1.0 else 0.0)
        T[i, i] <- c 
        T[j, j] <- c
        T[i, j] <- s 
        T[j, i] <- -s 
        let A_new = (Matrix.transpose T) * A * T
        if numIteration = maxIterations || nondiagonalSquaresSum A < epsilon then
            A_new, numIteration
        else
            step A_new (numIteration + 1)

    let (diagonalizedMatrix, numIterations) = step B 1
    let eigenValues = [0..n - 1] 
                      |> List.map (fun i -> diagonalizedMatrix[i, i])
    eigenValues, numIterations
            

let A1 = createHilbertMatrix 5
let eigenValues1 = [0.0; 0.0; 0.011; 0.209; 1.567]
let A2 = matrixFromList [[13.0; 16.0; 16.0];
                         [-5.0; -7.0; -6.0];
                         [-6.0; -8.0; -7.0]]
let eigenValues2 = [1.0; 1.0; -3.0]
let A3 = createHilbertMatrix 8
let eigenValues3 = [0.0; 0.0; 0.0; 0.0; 0.01; 0.026; 0.298; 1.696]
let A4 = matrixFromList [[1.0; 0.42; 0.54; 0.56]
                         [0.42; 1.00; 0.32; 0.44]
                         [0.54; 0.32; 1.0; 0.22]
                         [0.66; 0.44; 0.22; 1.0]]
let eigenValues4 = [0.242261; 0.638284; 0.796707; 2.322749]
let A5 = matrixFromList [[-5.0; 0.0; 0.0]
                         [0.0; -11.0; 0.0]
                         [0.0; 0.0; -12.0]]
let eigenValues5 = [-12.0; -11.0; -5.0]


let experiment () =
    let A = A5
    let eigenValues = eigenValues5
    let epsilonValues = [1e-1; 1e-2; 1e-3; 1e-4; 1e-5; 1e-6]
    let mutable (approximativeEigenValuesAbs, numIterationsAbs) = 
        JacobiMethod A chooseMaxAbsElement 1e-3 1000
    let mutable (approximativeEigenValuesGershgorin, numIterationsGershgorin) =
        JacobiMethod A chooseByMaxGershgorinCircle 1e-3 1000
    printfn "Приближённые значения: %A" (approximativeEigenValuesAbs |> List.sort)
    printfn "Приближённые значения: %A" (approximativeEigenValuesGershgorin |> List.sort)
    printfn "Точные значения: %A" eigenValues
    let check = GershgorinCirclesCheck A approximativeEigenValuesAbs &&
                GershgorinCirclesCheck A approximativeEigenValuesGershgorin
    match check with
    | true -> printfn "Собственные значения лежат в кругах Гершгорина"
    | false -> printfn "Собственные значения не лежат в кругах Гершгорина"

    printfn "epsilon  | diff. max nondiag. element | diff. max Gershgorin circle"

    let resultAbs = epsilonValues 
                    |> List.map (fun eps -> JacobiMethod A chooseMaxAbsElement eps 1000)
    let resultGershgorin = epsilonValues 
                           |> List.map (fun eps -> JacobiMethod A chooseByMaxGershgorinCircle eps 1000)
    
    for i in 0..epsilonValues.Length - 1 do
        let eigenValuesAbs = resultAbs[i] 
                             |> fst 
                             |> List.sort 
                             |> Array.ofList 
                             |> Vector.ofArray
        let numIterationsAbs = resultAbs[i] |> snd
        let eigenValuesGershgorin = resultGershgorin[i] 
                                    |> fst 
                                    |> List.sort 
                                    |> Array.ofList 
                                    |> Vector.ofArray
        let numIterationsGershgorin = resultGershgorin[i] |> snd
        let trueEigenValues = eigenValues |> Array.ofList |> Vector.ofArray
        let absDiff = eigenValuesAbs - trueEigenValues
                      |> Vector.norm
        let GershgorinDiff = eigenValuesGershgorin - trueEigenValues
                             |> Vector.norm
        printfn "%f   | %f     | %f" epsilonValues[i] absDiff GershgorinDiff
    done       
        

    let graphA = Chart.Point([1; 2; 3; 4; 5; 6], resultAbs |> List.map snd) 
                 |> Chart.withTitle "Max abs value strategy"
                 |> Chart.withXAxisStyle "epsilon"
                 |> Chart.withYAxisStyle "iterations"
    let graphB = Chart.Point([1; 2; 3; 4; 5; 6], resultGershgorin |> List.map snd) 
                 |> Chart.withTitle "Max Gershgorin circle strategy"
                 |> Chart.withXAxisStyle "epsilon"
                 |> Chart.withYAxisStyle "iterations"
    let combinedGraph = [graphA; graphB] |> Chart.combine
    combinedGraph |> Chart.show

experiment ()

