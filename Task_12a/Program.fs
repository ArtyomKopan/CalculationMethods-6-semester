open FsAlg.Generic


let mutable countFunctionCalls = 0

let f (x: float) (y: float) (z: float) =
    countFunctionCalls <- countFunctionCalls + 1
    (exp (2.0 * x)) * (sin (y * z)) * z

let f'x (x: float) (y: float) (z: float) =
    countFunctionCalls <- countFunctionCalls + 1
    2.0 * (exp (2.0 * x)) * (sin (y * z)) * z

let f'y (x: float) (y: float) (z: float) = 
    countFunctionCalls <- countFunctionCalls + 1
    (exp (2.0 * x)) * (cos (y * z)) * z * z

let f'z (x: float) (y: float) (z: float) =
    countFunctionCalls <- countFunctionCalls + 1
    (exp (2.0 * x)) * ((sin (y * z)) + (cos (y * z) * y * z))

let f''xx (x: float) (y: float) (z: float) =
    countFunctionCalls <- countFunctionCalls + 1
    4.0 * (exp (2.0 * x)) * (sin (y * z)) * z

let f''xy (x: float) (y: float) (z: float) =
    countFunctionCalls <- countFunctionCalls + 1
    2.0 * (exp (2.0 * x)) * (cos (y * z)) * z * z

let f''xz (x: float) (y: float) (z: float) =
    countFunctionCalls <- countFunctionCalls + 1
    2.0 * (exp (2.0 * x)) * ((sin (y * z)) + (cos (y * z) * y * z))

let f''yy (x: float) (y: float) (z: float) =
    countFunctionCalls <- countFunctionCalls + 1
    (exp (2.0 * x)) * -(sin (y * z)) * z * z * z

let f''yz (x: float) (y: float) (z: float) =
    countFunctionCalls <- countFunctionCalls + 1
    (exp (2.0 * x)) * (-y * z * z * (sin (y * z)) + 2.0 * z * (cos (y * z)))

let f''zz (x: float) (y: float) (z: float) = 
    countFunctionCalls <- countFunctionCalls + 1
    (exp (2.0 * x)) * (y * (cos (y * z)) + y * (cos (y * z)) - y * y * z * (sin (y * z)))

let secondDerivativesMatrix (x: float) (y: float) (z: float) =
    let J = Matrix.init 3 3 (fun i j -> 0.0)
    J[0, 0] <- f''xx x y z
    J[0, 1] <- f''xy x y z
    J[0, 2] <- f''xz x y z
    J[1, 0] <- f''xy x y z
    J[1, 1] <- f''yy x y z
    J[1, 2] <- f''yz x y z
    J[2, 0] <- f''xz x y z
    J[2, 1] <- f''yz x y z
    J[2, 2] <- f''zz x y z
    J

let gradient x y z =
    [f'x x y z; f'y x y z; f'z x y z] |> List.toArray |> Vector.ofArray

let gradientDescent epsilon maxIterations (alpha: float) =
    let stopCriterion (x: Vector<float>) =
        Vector.norm (gradient x[0] x[1] x[2]) < epsilon
        
    let rec step x num_iter =
        if num_iter = maxIterations || stopCriterion x then
            x, num_iter
        else
            let x_new = x - alpha * (gradient x[0] x[1] x[2]) 
            step x_new (num_iter + 1)

    let x0 = Vector.init 3 (fun i -> 1.0)
    step x0 1

let PolyakMethod epsilon maxIterations (alpha: float) (beta: float) =
    let stopCriterion (x: Vector<float>) =
        Vector.norm (gradient x[0] x[1] x[2]) < epsilon
        
    let rec step (x_prev: Vector<float>) (x: Vector<float>) num_iter =
        if num_iter = maxIterations || stopCriterion x then
            x, num_iter
        else
            let x_new = x - alpha * (gradient x[0] x[1] x[2]) + beta * (x - x_prev)
            step x x_new (num_iter + 1)

    let x0 = Vector.init 3 (fun i -> 1.0)
    let x1 = Vector.init 3 (fun i -> 2.0)
    step x0 x1 1

let NesterovMethod epsilon maxIterations (alpha: float) (beta: float) =
    let stopCriterion (x: Vector<float>) =
        Vector.norm (gradient x[0] x[1] x[2]) < epsilon

    let rec step (x: Vector<float>) (y: Vector<float>) num_iter =
        let x_new = y - alpha * (gradient y[0] y[1] y[2])
        if num_iter = maxIterations || stopCriterion x_new then
            x_new, num_iter
        else
            let y_new = x_new + beta * (x_new - x)
            step x_new y_new (num_iter + 1)

    let x0 = Vector.init 3 (fun i -> 1.0)
    let y0 = Vector.copy x0
    step x0 y0 1

let NewtonMethod epsilon maxIterations =
    let stopCriterion (x: Vector<float>) =
        Vector.norm (gradient x[0] x[1] x[2]) < epsilon
        
    let rec step x num_iter =
        if num_iter = maxIterations || stopCriterion x then
            x, num_iter
        else
            let HInv = secondDerivativesMatrix x[0] x[1] x[2] |> Matrix.inverse
            let x_new = x - HInv * (gradient x[0] x[1] x[2])  
            step x_new (num_iter + 1)

    let x0 = Vector.init 3 (fun i -> 1.0)
    step x0 1


let experiment () =
    let epsilon = 1e-6
    let maxIterations = 1000
    let alpha = 1.0
    let beta = 1.0

    let mutable totalTimeGradient = 0.0
    let mutable totalTimePolyak = 0.0
    let mutable totalTimeNesterov = 0.0
    let mutable totalTimeNewton = 0.0

    let mutable countGlobalMinsGradient = 0
    let mutable countGlobalMinsPolyak = 0
    let mutable countGlobalMinsNesterov = 0
    let mutable countGlobalMinsNewton = 0

    let mutable countFunctionCallsGradient = 0
    let mutable countFunctionCallsPolyak = 0
    let mutable countFunctionCallsNesterov = 0
    let mutable countFunctionCallsNewton = 0

    let mutable accuracyGradient = 0.0
    let mutable accuracyPolyak = 0.0
    let mutable accuracyNesterov = 0.0
    let mutable accuracyNewton = 0.0

    for i in 1..10 do
        let timer = System.Diagnostics.Stopwatch.StartNew()
        let x, numIters = gradientDescent (epsilon / 10.0 ** i) maxIterations alpha
        timer.Stop()
        totalTimeGradient <- totalTimeGradient + float(timer.ElapsedMilliseconds)
        accuracyGradient <- accuracyGradient + abs ((f x[0] x[1] x[2]) - 0.0)
        if abs ((f x[0] x[1] x[2]) - 0.0) < epsilon / 10.0 ** i then 
            countGlobalMinsGradient <- countGlobalMinsGradient + 1
        countFunctionCallsGradient <- countFunctionCallsGradient + countFunctionCalls
        countFunctionCalls <- 0
    done

    for i in 1..10 do
        let timer = System.Diagnostics.Stopwatch.StartNew()
        let x, numIters = PolyakMethod (epsilon / 10.0 ** i) maxIterations alpha beta
        timer.Stop()
        totalTimePolyak <- totalTimePolyak + float(timer.ElapsedMilliseconds)
        accuracyPolyak <- accuracyPolyak + abs ((f x[0] x[1] x[2]) - 0.0)
        if abs ((f x[0] x[1] x[2]) - 0.0) < epsilon / 10.0 ** i then 
            countGlobalMinsPolyak <- countGlobalMinsPolyak + 1
        countFunctionCallsPolyak <- countFunctionCallsPolyak + countFunctionCalls
        countFunctionCalls <- 0
    done

    for i in 1..10 do
        let timer = System.Diagnostics.Stopwatch.StartNew()
        let x, numIters = NesterovMethod (epsilon / 10.0 ** i) maxIterations alpha beta
        timer.Stop()
        totalTimeNesterov <- totalTimeNesterov + float(timer.ElapsedMilliseconds)
        accuracyNesterov <- accuracyNesterov + abs ((f x[0] x[1] x[2]) - 0.0)
        if abs ((f x[0] x[1] x[2]) - 0.0) < epsilon / 10.0 ** i then 
            countGlobalMinsNesterov <- countGlobalMinsNesterov + 1
        countFunctionCallsNesterov <- countFunctionCallsNesterov + countFunctionCalls
        countFunctionCalls <- 0
    done

    for i in 1..10 do
        let timer = System.Diagnostics.Stopwatch.StartNew()
        let x, numIters = gradientDescent (epsilon / 10.0 ** i) maxIterations alpha
        timer.Stop()
        totalTimeNewton <- totalTimeNewton + float(timer.ElapsedMilliseconds)
        accuracyNewton <- accuracyNewton + abs ((f x[0] x[1] x[2]) - 0.0)
        if abs ((f x[0] x[1] x[2]) - 0.0) < epsilon / 10.0 ** i then 
            countGlobalMinsNewton <- countGlobalMinsNewton + 1
        countFunctionCallsNewton <- countFunctionCallsNewton + countFunctionCalls
        countFunctionCalls <- 0
    done

    let meanTimeGradient = System.Math.Round(totalTimeGradient / 10.0, 5)
    let meanTimePolyak = System.Math.Round(totalTimePolyak / 10.0, 5)
    let meanTimeNesterov = System.Math.Round(totalTimeNesterov / 10.0, 5)
    let meanTimeNewton = System.Math.Round(totalTimeNewton / 10.0, 5)

    let meanAccuracyGradient = System.Math.Round(accuracyGradient / 10.0, 5)
    let meanAccuracyPolyak = System.Math.Round(accuracyPolyak / 10.0, 5)
    let meanAccuracyNesterov = System.Math.Round(accuracyNesterov / 10.0, 5)
    let meanAccuracyNewton = System.Math.Round(accuracyNesterov / 10.0, 5)

    printfn "Method  |  Time  |  Calls  |  Accuracy  |  Global mins\n"
    printfn "Gradient| %f  | %d  | %f  | %d\n" meanTimeGradient countFunctionCallsGradient meanAccuracyGradient countGlobalMinsGradient
    printfn "Polyak  | %f  | %d  | %f  | %d\n" meanTimePolyak countFunctionCallsPolyak meanAccuracyPolyak countGlobalMinsPolyak
    printfn "Nesterov| %f  | %d  | %f  | %d\n" meanTimeNesterov countFunctionCallsNesterov meanAccuracyNesterov countGlobalMinsNesterov
    printfn "Newton  | %f  | %d  | %f  | %d\n" meanTimeNewton countFunctionCallsNewton meanAccuracyNewton countGlobalMinsNewton

experiment ()