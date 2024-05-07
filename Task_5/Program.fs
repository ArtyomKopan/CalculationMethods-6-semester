open FsAlg.Generic
open Plotly.NET

let createHilbertMatrix n =
    Matrix.init n n (fun i j -> 1.0 / (float(i) + float(j) + 1.0))

let matrixFromList inputList =
    let numRows = List.length inputList
    let numCols = List.length (List.head inputList)
    let array2D = Array2D.init numRows numCols (fun i j -> inputList.[i].[j])
    Matrix.ofArray2D array2D

let dotProduct (x: Vector<float>) (y: Vector<float>) =
    List.sumBy (fun i -> x[i] * y[i]) [0..x.Length - 1]

let normalize (x: Vector<float>) =
    let norm = Vector.norm x
    if norm > 1.0 then
        (1.0 / norm) * x
    else
        x

let powerMethod (A: Matrix<float>) epsilon maxIterations =
    let n = A.Rows
    let x0 = Vector.init n (fun i -> 1.0)
    let eigenValue0 = 1.0
    
    let rec step currentEigenValue (x: Vector<float>) numIteration =
        let x_new = A * x
        let newEigenValue = (dotProduct x_new x_new) / (dotProduct x x) |> sqrt
        if numIteration = maxIterations || abs (newEigenValue - currentEigenValue) < epsilon then
            newEigenValue, normalize x_new, numIteration
        else
            step newEigenValue (normalize x_new) (numIteration + 1)
    
    step eigenValue0 x0 1

let dotProductsMethod (A: Matrix<float>) epsilon maxIterations =
    let n = A.Rows
    let ATransposed = A |> Matrix.transpose
    let x0 = Vector.init n (fun i -> 1.0)
    let y0 = Vector.copy x0
    let eigenValue0 = 1.0

    let rec step currentEigenValue (x: Vector<float>) (y: Vector<float>) numIteration =
        let x_new = A * x |> normalize
        let y_new = ATransposed * y |> normalize
        let newEigenValue = (dotProduct (A * x) y_new) / (dotProduct x y_new)
        if numIteration = maxIterations || abs (newEigenValue - currentEigenValue) < epsilon then
            newEigenValue, x_new, y_new, numIteration
        else
            step newEigenValue x_new y_new (numIteration + 1)

    step eigenValue0 x0 y0 1

let minPowerMethod (A: Matrix<float>) epsilon maxIterations =
    powerMethod (-A) epsilon maxIterations

let minDotProductsMethod (A: Matrix<float>) epsilon maxIterations =
    dotProductsMethod (-A) epsilon maxIterations


let A1 = createHilbertMatrix 5
let maxEigenValue1 = 1.567051
let A2 = matrixFromList [[3.27; 1.04; -1.37];
                         [1.04; 2.97; 0.93];
                         [1.37; 0.93; 4.83]]
let maxEigenValue2 = 4.5227
let A3 = createHilbertMatrix 10
let maxEigenValue3 = 1.751920
let A4 = matrixFromList [[-5.51; 1.87; 0.42]
                         [0.29; -11.81; 5.71]
                         [0.05; 4.31; -12.97]]
let maxEigenValue4 = 17.39
let A5 = matrixFromList [[-5.0; 0.0; 0.0]
                         [0.0; -11.0; 0.0]
                         [0.0; 0.0; -12.0]]
let maxEigenValue5 = 5.0


let experiment () =
    let A = A5
    let maxEigenValue = maxEigenValue5
    let epsilonValues = [1e-1; 1e-2; 1e-3; 1e-4; 1e-5; 1e-6]
    let resultPower = List.map (fun eps -> 
                                let (_, _, numIters) = powerMethod A eps 1000
                                numIters) epsilonValues
    let resultDotProducts = List.map (fun eps -> 
                                let (_, _, _, numIters) = dotProductsMethod A eps 1000
                                numIters) epsilonValues

    let powerEigenValues = List.map (fun eps -> 
                                       let (eigenValue, _, _) = powerMethod -A eps 1000
                                       eigenValue) epsilonValues
    let dotProductsEigenValues = List.map (fun eps -> 
                                       let (eigenValue, _, _, _) = dotProductsMethod -A eps 1000
                                       eigenValue) epsilonValues
    //let maxEigenValue = Matrix.eigenvalues A |> Vector.toArray |> Array.max
    printfn "%f\n" maxEigenValue
    printfn "%f\n" powerEigenValues[5]
    printfn "%f\n" dotProductsEigenValues[5]
    printfn "epsilon    | diff. power method | diff. dot products method"
    for i in 0..epsilonValues.Length - 1 do
        let powerDifference = abs (powerEigenValues[i] - maxEigenValue)
        let dotProductsDifference = abs (dotProductsEigenValues[i] - maxEigenValue)
        printfn "%f     | %f                 | %f" epsilonValues[i] powerDifference dotProductsDifference
    done

    let graphA = Chart.Point([1; 2; 3; 4; 5; 6], resultPower) 
                 |> Chart.withTitle "Power method"
                 |> Chart.withXAxisStyle "epsilon"
                 |> Chart.withYAxisStyle "iterations"
    let graphB = Chart.Point([1; 2; 3; 4; 5; 6], resultDotProducts) 
                 |> Chart.withTitle "Dot products method"
                 |> Chart.withXAxisStyle "epsilon"
                 |> Chart.withYAxisStyle "iterations"
    let combinedGraph = [graphA; graphB] |> Chart.combine
    combinedGraph |> Chart.show

experiment ()

