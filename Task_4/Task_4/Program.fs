open FsAlg.Generic

let rnd = new System.Random()

let matrixNorm (a: Matrix<float>) N = 
    let sumOfAbs j =
        List.fold (fun acc i -> acc + abs a[i, j]) 0.0 [0..N-1]
    List.init N sumOfAbs |> List.max

let createDiagonalMatrix n =
    Matrix.init n n (fun i j -> if i = j then rnd.NextDouble() else 0.0)

let solveEqutationWithDiagonalMatrix (A: Matrix<float>) (b: Vector<float>) =
    let n = A.Rows
    List.map (fun i -> b[i] / A[i, i]) [0..n - 1] 
    |> List.toArray 
    |> Vector.ofArray
    
let matrixEquals (A: Matrix<float>) (B: Matrix<float>) epsilon =
    let mutable nonEqual = false
    if A.Rows <> B.Rows || A.Cols <> B.Cols then
        nonEqual <- true
    if nonEqual then false
    else 
        for i in 0..A.Rows - 1 do
            for j in 0..A.Cols - 1 do 
                if abs (A[i, j] - B[i, j]) > epsilon || nonEqual then
                    nonEqual <- true
                else
                    nonEqual <- false
            done
        done
        not nonEqual

let isDiagonalPredominance (A: Matrix<float>) =
    let n = A.Rows
    
    let rowDiagonalPredominance i =
        let nonDiagonalSum = (Matrix.row i A 
                            |> Matrix.toVector 
                            |> Vector.toArray 
                            |> Array.sum) - A[i, i]
        abs nonDiagonalSum < abs A[i, i]
        
    let columnDiagonalPredominance i =
        let nonDiagonalSum = (Matrix.col i A
                            |> Matrix.toVector
                            |> Vector.toArray
                            |> Array.sum) - A[i, i]
        abs nonDiagonalSum < abs A[i, i]
    
    [0..n - 1] 
    |> List.map (fun i -> rowDiagonalPredominance i && columnDiagonalPredominance i) 
    |> List.forall (fun item -> item)

let createRandomVector n =
    Vector.init n (fun i -> rnd.NextDouble())
    
let createDiagonalDominanceMatrix dimension =
    let diagonalElements = Array.init dimension (fun _ -> rnd.NextDouble() * 10000.0)
    Matrix.init dimension dimension
        (fun i j ->
        if i = j
        then diagonalElements[i]
        else diagonalElements[i] / float(dimension + 1))

let areSolutionsNear (x: Vector<float>) (y: Vector<float>) eps =
    Vector.norm (x - y) < eps

let Seidel (A: Matrix<float>) (b: Vector<float>) epsilon maxIterations =
    let n = A.Rows
    let D = Matrix.init n n (fun i j -> if i = j then A[i, i] else 0.0)
    let L = Matrix.init n n (fun i j -> if i > j then A[i, j] else 0.0)
    let R = Matrix.init n n (fun i j -> if i < j then A[i, j] else 0.0)
    let multiplierMatrix = D + L |> Matrix.inverse
    let A_norm = matrixNorm A n
    
    let rec iterate (x_prev: Vector<float>) k =
        let x_new = -multiplierMatrix * R * x_prev + multiplierMatrix * b
        let diff = (A_norm / (1.0 - A_norm)) * Vector.norm (x_new - x_prev)
        if diff < epsilon || k = maxIterations then
            x_new, k
        else
            iterate x_new (k + 1)
    
    iterate (Vector.init n (fun i -> 0.0)) 1

let calcIterations (B: Matrix<float>) (c: Vector<float>) epsilon maxIterations =
    let n = B.Rows
    let normB = matrixNorm B n
    let rec iter (predX: Vector<float>) (currX: Vector<float>) numIter =
        let newX = B * currX + c
        let diff = normB * Vector.norm (newX - currX) / (1.0 - normB)
        if numIter = maxIterations || abs diff < epsilon then
            newX, numIter
        else
            iter currX newX (numIter + 1)
    iter (Vector.init n (fun i -> 0.0)) (Vector.copy c) 0

let simpleIterationMethod (A: Matrix<float>) (b: Vector<float>) epsilon maxIterations =
    let n = A.Rows        
    let B = Matrix.init n n (fun i j -> if i <> j then -A[i, j] / A[i, i] else 0.0)
    let c = Vector.init n (fun i -> b[i] / A[i, i])
    if (Matrix.eigenvalues B |> Vector.toArray |> Array.max) > matrixNorm B n then
        failwith "Max eigen value > ||B||"
    calcIterations B c epsilon maxIterations

let runTests tests =
    let runTest test = 
        let isPassed = test()
        match isPassed with
        | true -> 
            printfn "TEST PASSED"
            isPassed
        | false ->
            printfn "TEST FAILED"
            isPassed
    List.forall (fun test -> runTest test) tests

let SeidelTest1 () =
    printfn "Тест 1. Метод Зейделя"
    let A = createDiagonalMatrix 200
    let b = createRandomVector 200
    let eps = 1e-4
    let maxIterations = 1000
    let x0 = solveEqutationWithDiagonalMatrix A b 
    let (x, nIterations) = Seidel A b eps maxIterations
    printfn "Число итераций = %d" nIterations
    areSolutionsNear x x0 1e-3
    
let SeidelTest2 () = 
    printfn "Тест 2. Метод Зейделя"
    let A = createDiagonalDominanceMatrix 50
    let b = createRandomVector 50
    let eps = 1e-5
    let maxIterations = 1000
    let x0 = Matrix.solve A b
    let (x, nIterations) = Seidel A b eps maxIterations
    printfn "Число итераций = %d" nIterations
    areSolutionsNear x x0 1e-3

let simpleIterationMethodTest1 () = 
    printfn "Тест 3. Метод простой итерации"
    let A = createDiagonalMatrix 50
    let b = createRandomVector 50
    let eps = 1e-4
    let maxIterations = 1000
    let x0 = solveEqutationWithDiagonalMatrix A b 
    let (x, nIterations) = simpleIterationMethod A b eps maxIterations
    printfn "Число итераций = %d" nIterations
    areSolutionsNear x x0 1e-3
    
let simpleIterationMethodTest2 () = 
    printfn "Тест 4. Метод простой итерации"
    let A = createDiagonalDominanceMatrix 50
    let b = createRandomVector 50
    let eps = 1e-4
    let maxIterations = 1000
    let x0 = Matrix.solve A b
    let (x, nIterations) = simpleIterationMethod A b eps maxIterations
    printfn "Число итераций = %d" nIterations
    areSolutionsNear x x0 1e-3

runTests [SeidelTest1; SeidelTest2; simpleIterationMethodTest1; simpleIterationMethodTest2] |> ignore