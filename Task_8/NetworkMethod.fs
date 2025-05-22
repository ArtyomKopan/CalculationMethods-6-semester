module NetworkMethod

open FsAlg.Generic
open Plotly.NET


// граничные условия имеют вид f(x0) = y0
type Condition(x0: float, y0: float) =
    member c.x = x0
    member c.y = y0

// уравнение вида -p(x)u'' + q(x)u' + r(x)u = f(x)
type Equation(p: (float -> float),
              q: (float -> float),
              r: (float -> float),
              f: (float -> float), 
              condition1: Condition, 
              condition2: Condition) =
    member e.p = p
    member e.q = q
    member e.r = r
    member e.f = f
    member e.condition1 = condition1
    member e.condition2 = condition2

let buildNetwork (a: float) (b: float) (n: int) =
    let h = (b - a) / float(n)
    List.init (n + 1) (fun i -> float(a) + float(i) * h)

let solveWithNetwork (equation: Equation) a b n0 epsilon maxIterations =
    let alpha = equation.condition1.y
    let beta = equation.condition2.y
    let r = 2.0
    let p = 1.0
    let mutable errors = []
    let mutable networkSteps = [n0]

    let error (v1: Vector<float>) (v2: Vector<float>) =
        let delta = Vector.init 
                                v2.Length 
                                (fun n -> if n % 2 = 0 then (v2[n] - v1[n / 2]) / (r ** p - 1.0) else 0.0)
        for n in 1..2..v2.Length - 1 do
            delta[n] <- (delta[n - 1] + delta[n + 1]) / 2.0
        done
        delta

    let errorNorm (v1: Vector<float>) (v2: Vector<float>) =
        Vector.norm (error v1 v2) / (v2.Length |> float |> sqrt)

    // уточняем решение v2
    let clarificate (v1: Vector<float>) (v2: Vector<float>) = v2 + error v1 v2

    let rec step n_prev (v_prev: Vector<float>) numIteration =
        printfn "Iteration %d" numIteration
        let n = n_prev * 2
        let h = (b - a) / float(n)
        let x = buildNetwork a b n
        let p = List.init (n + 1) (fun i -> equation.p x[i])
        let q = List.init (n + 1) (fun i -> equation.q x[i])
        let r = List.init (n + 1) (fun i -> equation.r x[i])
        let f = Vector.init (n + 1) (fun i -> equation.f x[i])
        // нужно получить СЛАУ Az = f
        let A = Matrix.init (n + 1) (n + 1) (fun i j -> 0.0)
        A[0, 0] <- 1.0
        A[n, n] <- 1.0
        f[0] <- alpha
        f[n] <- beta
        
        for i in 1..n - 1 do
            A[i, i + 1] <- p[i] / (h ** 2.0) + q[i] / (2.0 * h)
            A[i, i] <- -2.0 * p[i] / (h ** 2.0) - r[i]
            A[i, i - 1] <- p[i] / (h ** 2.0) - q[i] / (2.0 * h)
        done

        let v_current = Matrix.solve A f
        let difference = errorNorm v_prev v_current
        errors <- difference :: errors
        networkSteps <- n :: networkSteps
        
        if difference < epsilon || numIteration = maxIterations then 
            clarificate v_prev v_current
        else 
            step n v_current (numIteration + 1)

    let v0 = Vector.init (n0 + 1) (fun i -> 0.0)
    let solution = step n0 v0 1
    solution, errors, networkSteps        

let drawPlot (x: List<float>) (trueSolution: float -> float) (approxSolution: Vector<float>) =
    if x.Length <> approxSolution.Length then
        failwith "Length of x vector must be equal length of approxSolution vector"
    let y = x |> List.map (fun xi -> trueSolution xi) |> List.toSeq
    let xy_true = Seq.zip x y
    let xy_approx = Seq.zip x (approxSolution |> Vector.toSeq)
    let truePlot = Chart.Line xy_true 
                   |> Chart.withTitle "Сравнение точного и приближённого решения"
                   |> Chart.withXAxisStyle "x"
                   |> Chart.withYAxisStyle "u(x)"
    let approxPlot = Chart.Line xy_approx
    let plot = [truePlot; approxPlot] |> Chart.combine
    plot |> Chart.show

let drawAccuracyPlot (n: List<int>) (accuracy: List<float>) =
    let xy = Seq.zip (n |> List.rev |> List.toSeq) (accuracy |> List.rev |> List.toSeq)
    let plot = Chart.Line xy
               |> Chart.withTitle "Зависимость погрешности от количества итераций"
               |> Chart.withXAxisStyle "Количество итераций"
               |> Chart.withYAxisStyle "Погрешность"
    plot |> Chart.show


// вариант 3 из методички
let test3 () =
    let p x = 1.0 / (x - 3.0)
    let q x = 1.0 + x / 2.0
    let r x = -exp (x / 2.0)
    let f x = 2.0 - x
    let condition1 = Condition(-1.0, 0.0)
    let condition2 = Condition(1.0, 0.0)
    let equation = Equation(p, q, r, f, condition1, condition2)
    
    let n0 = 10
    let a = -1.0
    let b = 1.0
    let (solution, errors, networkSteps) = solveWithNetwork equation a b n0 1e-6 6
    let n = solution.Length
    let x = buildNetwork a b (n - 1)
    let u x = 0.0 // Fake solution
    drawPlot x u solution
    drawAccuracyPlot networkSteps errors

let test1 () =
    // u'' = 1, u(0) = 0, u(1) = 1
    // u = x^2/2 + 1/2*x
    let p x = 1.0
    let q x = 0.0
    let r x = 0.0
    let f x = 1.0
    let condition1 = Condition(0.0, 0.0)
    let condition2 = Condition(1.0, 1.0)
    let equation = Equation(p, q, r, f, condition1, condition2)
    let u x = 0.5 * x ** 2.0 + 0.5 * x

    let n0 = 10
    let a = 0.0
    let b = 1.0
    let (solution, errors, networkSteps) = solveWithNetwork equation a b n0 1e-6 7
    let n = solution.Length
    let x = buildNetwork a b (n - 1)
    let expected = Vector.init n (fun i -> u x[i])
    printfn "%A" <| Vector.toArray solution
    printfn "%f" <| Vector.norm (expected - solution)
    drawPlot x u solution
    drawAccuracyPlot networkSteps errors
    
let test2 () =
    // (x-1)u'' - xu' + u = (x - 1)^2
    // u(x) = -1 + exp(x) + 5x/2 - 2sqrt(e)x - x^2
    let p x = x - 1.0
    let q x = -x
    let r x = 1.0
    let f x = (x - 1.0) ** 2.0
    let u x = -1.0 + exp x + 2.5 * x - 2.0 * (sqrt System.Math.E) * x - x ** 2.0
    let condition1 = Condition(0.0, 0.0)
    let condition2 = Condition(0.5, 0.0)
    let equation = Equation(p, q, r, f, condition1, condition2)

    let n0 = 10
    let a = 0.0
    let b = 0.5
    let (solution, errors, networkSteps) = solveWithNetwork equation a b n0 1e-6 7
    let n = solution.Length
    let x = buildNetwork a b (n - 1)
    let expected = Vector.init n (fun i -> u x[i])
    printfn "%A" <| Vector.toArray solution
    printfn "%f" <| Vector.norm (expected - solution)
    drawPlot x u solution
    drawAccuracyPlot networkSteps errors

let test4 () =
    // u'' - 2u' + 10u = 0
    // u(x) = exp(x - 1) / sin(6)  * (14 * sin(3x+3) - exp(2) * sin(3 - 3x))
    let p x = 1.0
    let q x = -2.0
    let r x = 10.0
    let f x = 0.0
    let u x = (exp (x - 1.0)) / (sin 6.0) * (14.0 * sin (3.0 * x + 3.0) - (exp 2.0) * sin (3.0 - 3.0 * x))
    let condition1 = Condition(-1.0, -1.0)
    let condition2 = Condition(1.0, 14.0)
    let equation = Equation(p, q, r, f, condition1, condition2)

    let n0 = 10
    let a = -1.0
    let b = 1.0
    let (solution, errors, networkSteps) = solveWithNetwork equation a b n0 1e-6 7
    let n = solution.Length
    let x = buildNetwork a b (n - 1)
    let expected = Vector.init n (fun i -> u x[i])
    printfn "%A" <| Vector.toArray solution
    printfn "%f" <| Vector.norm (expected - solution)
    drawPlot x u solution
    drawAccuracyPlot networkSteps errors

let test5 () =
    // u'' - 2u' + y = 0
    // u(x) = 0.5 * exp(x - 1) * (e ** 2 * (x - 1) + 14 * (x + 1))
    let p x = 1.0
    let q x = -2.0
    let r x = 1.0
    let f x = 0.0
    let u x = 0.5 * exp (x - 1.0) * (System.Math.E ** 2.0 * (x - 1.0) + 14.0 * (x + 1.0))
    let condition1 = Condition(-1.0, -1.0)
    let condition2 = Condition(1.0, 14.0)
    let equation = Equation(p, q, r, f, condition1, condition2)

    let n0 = 10
    let a = -1.0
    let b = 1.0
    let (solution, errors, networkSteps) = solveWithNetwork equation a b n0 1e-6 7
    let n = solution.Length
    let x = buildNetwork a b (n - 1)
    let expected = Vector.init n (fun i -> u x[i])
    printfn "%A" <| Vector.toArray solution
    printfn "%f" <| Vector.norm (expected - solution)
    drawPlot x u solution
    drawAccuracyPlot networkSteps errors