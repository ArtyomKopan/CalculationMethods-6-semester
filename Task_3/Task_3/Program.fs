open FsAlg.Generic


let createHilbertMatrix n =
    Matrix.init n n (fun i j -> 1.0 / (float(i) + float(j) + 1.0))

let cond_v (A: Matrix<double>) =
    let N = A.Rows
    let squares_sum n =
        List.fold (fun acc m -> acc + A.[n, m] ** 2) 0.0 [0..N-1]
        |> sqrt
    let product = List.fold (fun acc n -> acc * (squares_sum n)) 1.0 [0..N-1]
    product / (Matrix.det A)

let cond_a (A: Matrix<double>) =
    let N = A.Rows
    let C = Matrix.inverse A
    List.map 
        (fun i -> Vector.norm (Matrix.toVector A[i,*]) * Vector.norm (Matrix.toVector C[*,i])) 
        [0..N-1]
    |> List.max

let matrixNorm (a: Matrix<double>) (N: int) =
    let sumOfAbs j =
        List.fold (fun acc i -> acc + abs a.[i, j]) 0.0 [0..N-1]
    List.init N sumOfAbs 
    |> List.max

let cond_s (A: Matrix<double>) =
    let N = A.Rows
    (matrixNorm A N) * (matrixNorm (Matrix.inverse A) N)

let rotationMethod (A: Matrix<float>) =
    let n = A.Rows
    let mutable Q = Matrix.init n n (fun i j -> if i = j then 1.0 else 0.0)
    let mutable R = Matrix.copy A
    
    for i in 0..n - 1 do
        for j in i + 1..n - 1 do
            let T = Matrix.init n n (fun i j -> if i = j then 1.0 else 0.0)
            let square_root = R[i, i] ** 2 + R[j, i] ** 2 |> sqrt
            let c = R[i, i] / square_root
            let s = -R[j, i] / square_root
            T[i, i] <- c
            T[j, i] <- s
            T[j, j] <- c
            T[i, j] <- -s
            R <- T * R
            Q <- T * Q
        done
    done
    
    for i in 0..n - 1 do
        for j in i + 1..n - 1 do
            R[j, i] <- 0
        done
    done
    
    Matrix.inverse Q, R

let solveWithSubstitution (A: Matrix<float>) (b: Vector<float>) =
    let n = A.Rows
    let x = Vector.init n (fun i -> 0.0)
    for k in 0..n - 1 do
        x[k] <- (-) b[k] <| List.sumBy (fun i -> A[k, i] * x[i]) [0..k - 1]
        x[k] <- x[k] / A[k, k]
    done
    x

let solveWithQR (A: Matrix<float>) (b: Vector<float>) =
    let Q, R = rotationMethod A
    let y = (Matrix.transpose Q) * b
    let x = Matrix.solve R y
    x

type Test(id: int, name: string, test_function: (unit -> bool)) =
    member test.id = id
    member test.name = name
    member test.test_function = test_function

let runTests (tests: List<Test>) =
    let runTest (test: Test)=
        let isPassed = test.test_function()
        if isPassed then
            printfn "TEST %d PASSED (%s)" test.id test.name
        else
            printfn "TEST %d FAILED (%s)" test.id test.name
        isPassed
    
    List.forall (fun test -> runTest test) tests

let testEquation1 () =
    let A = createHilbertMatrix 5
    let x = Vector.init 5 (fun i -> 1.0)
    let b = A * x
    let solution = solveWithQR A b
    Vector.norm (x - solution) < 0.1
    
let testEquation2 () =
    let A = createHilbertMatrix 8
    let x = Vector.init 8 (fun i -> 1.0)
    let b = A * x
    let solution = solveWithQR A b
    Vector.norm (x - solution) < 0.1
    
let testEquation3 () =
    let inputList = [[3.27; 1.04; -1.37];
                 [1.04; 2.97; 0.93];
                 [1.37; 0.93; 4.83]]

    let numRows = List.length inputList
    let numCols = List.length (List.head inputList)

    let array2D = Array2D.init numRows numCols (fun i j -> inputList.[i].[j])

    let A = Matrix.ofArray2D array2D
    let x = Vector.init 3 (fun i -> 1.0)
    let b = A * x
    let solution = solveWithQR A b
    Vector.norm (x - solution) < 0.1
    
let testConds1 () =
    let A = createHilbertMatrix 5
    let Q, R = rotationMethod A
    printfn "cond_v A = %f, cond_a A = %f, cond_s A = %f\n" (cond_v A) (cond_a A) (cond_s A)
    printfn "cond_v Q = %f, cond_a Q = %f, cond_s Q = %f\n" (cond_v Q) (cond_a Q) (cond_s Q)
    printfn "cond_v R = %f, cond_a R = %f, cond_s R = %f\n" (cond_v R) (cond_a R) (cond_s R)
    
let testConds2 () =
    let inputList = [[3.27; 1.04; -1.37];
                 [1.04; 2.97; 0.93];
                 [1.37; 0.93; 4.83]]

    let numRows = List.length inputList
    let numCols = List.length (List.head inputList)

    let array2D = Array2D.init numRows numCols (fun i j -> inputList.[i].[j])

    let A = Matrix.ofArray2D array2D
    let Q, R = rotationMethod A
    printfn "cond_v A = %f, cond_a A = %f, cond_s A = %f\n" (cond_v A) (cond_a A) (cond_s A)
    printfn "cond_v Q = %f, cond_a Q = %f, cond_s Q = %f\n" (cond_v Q) (cond_a Q) (cond_s Q)
    printfn "cond_v R = %f, cond_a R = %f, cond_s R = %f\n" (cond_v R) (cond_a R) (cond_s R)
