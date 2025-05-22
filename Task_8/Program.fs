open FsAlg.Generic
open Plotly.NET
open FSharp.Stats.SpecialFunctions.Gamma
open FSharp.Stats.Integration.Differentiation.TwoPointDifferentiation
open NetworkMethod
open System
open MathNet.Numerics.Differentiation
open MathNet.Symbolics

let integrate (f: float -> float) a b mode =
    let left_rectangle () =
        let x = buildNetwork a b 100
        List.sumBy (fun i -> (f x[i]) * (x[i + 1] - x[i])) [0..99]
        
    let right_rectangle () = (f b) * (b - a)

    let middle_rectangle () = (f (a + b) / 2.0) * (b - a)

    let trapezoids () = (f a + f b) / 2.0 * (b - a)

    let parabola () = ((b - a) / 6.0) * (f a + 4.0 * f ((a + b) / 2.0) + f b)
        
    match mode with 
    | "left rectangle" -> left_rectangle ()
    | "right rectangle" -> right_rectangle ()
    | "middle rectangle" -> middle_rectangle ()
    | "trapezoids" -> trapezoids ()
    | "parabola" -> parabola ()
    | _ -> left_rectangle ()

// скалярное произведение функций на отрезке [a, b]: (f, g) = интеграл от a до b (fg)dx
let dot (f: float -> float) (g: float -> float) a b mode = 
    integrate (fun x -> f x * g x) a b mode

let rec factorial n =
    if n <= 1 then 1
    else n * factorial (n - 1)

// многочлен Якоби порядка n при alpha = beta = 1
let JacobiPolynom n alpha beta x =
    // binomial coefficient
    let bc (z: float) (m: int) =
        if m < 0 then 0.0
        else gamma (z + 1.0) / (gamma (float(m) + 1.0) * gamma (z - float(m) + 1.0))

    let term s =
        (bc (float (n) + alpha) (n - s)) * 
        (bc (float (n) + beta) s) * 
        (((x - 1.0) / 2.0) ** float(s)) * 
        (((x + 1.0) / 2.0) ** float(n - s))

    List.sumBy (fun s -> term s) [0..n]

let w_i1 i x = (1.0 - x ** 2.0) ** 2.0 * (JacobiPolynom (i - 2) 2.0 2.0 x)

let w_i1_prime i x = -2.0 * float(i - 1) * (1.0 - x ** 2.0) * (JacobiPolynom (i - 1) 1.0 1.0 x)

let w_i1_second i x = 4.0 * float(i - 1) * float(i) * (JacobiPolynom i 0.0 0.0 x)

// k-ая производная многочлена Якоби
let JacobiDerivative n alpha beta k x =
    (gamma (alpha + beta + float(n) + 1.0 + float(k))) / 
    ((2.0 ** float(k)) * gamma (alpha + beta + float(n) + 1.0)) * 
    (JacobiPolynom (n - k) (alpha + float(k)) (beta + float(k)) x)

let GalerkinMethod (equation: Equation) a b n =
    let L i x = 
        (equation.p x) * (JacobiDerivative i 1.0 1.0 2 x) +    
        (equation.q x) * (JacobiDerivative i 1.0 1.0 1 x) +
        (equation.r x) * (JacobiPolynom i 1.0 1.0 x)

    //let L i x =     
    //    (equation.p x) * (w_i1_second i x) +
    //    (equation.q x) * (w_i1_prime i x) +
    //    (equation.r x) * (w_i1 i x)

    let w = List.map (fun i -> JacobiPolynom i 1.0 1.0) [0..n - 1]
    let Lw = [0..n - 1] |> List.map (fun i -> L i)

    // нужно получить систему уравнений Ac = b
    let A = Matrix.init n n (fun i j -> dot Lw[j] w[i] a b "left_rectangle")
    let b = Vector.init n (fun i -> dot equation.f w[i] a b "left_rectangle")
    for i in 0..n - 1 do
        printfn "%A" (Matrix.row i A)

    let c = Matrix.solve A b
    let u x = List.sumBy (fun i -> c[i] * (w[i] x)) [0..n - 1]
    u

let drawPlot3 (x: List<float>) (expected: List<float>) (networkSolution: List<float>) (GalerkinSolution: List<float>) =
    if x.Length <> expected.Length then
        failwith "Length of x != Length of expected"
    if x.Length <> networkSolution.Length then
        failwith "Length of x != Length of networkSolution"
    if x.Length <> GalerkinSolution.Length then
        failwith "Length of x != Length of GalerkinSolution"

    let xyTrue = Seq.zip x expected
    let xyNetwork = Seq.zip x networkSolution
    let xyGalerkin = Seq.zip x GalerkinSolution
    let truePlot = Chart.Line xyTrue
                   |> Chart.withTitle "Сравнение точного и приближённого решения"
                   |> Chart.withXAxisStyle "x"
                   |> Chart.withYAxisStyle "u(x)"
    let networkPlot = Chart.Line xyNetwork
    let GalerkinPlot = Chart.Line xyGalerkin
    let plot = [truePlot; networkPlot; GalerkinPlot] |> Chart.combine
    plot |> Chart.show

let test1 () =
    // u'' = 1, u(0) = 0, u(1) = 1
    // u = x^2/2 + 1/2*x
    let p x = 1.0
    let q x = -2.0
    let r x = 10.0
    let f x = 0.0
    let u x = (exp (x - 1.0)) / (sin 6.0) * (14.0 * sin (3.0 * x + 3.0) - (exp 2.0) * sin (3.0 - 3.0 * x))
    let condition1 = Condition(-1.0, -1.0)
    let condition2 = Condition(1.0, 14.0)
    let equation = Equation(p, q, r, f, condition1, condition2)

    let n = 10
    let a = -1.0
    let b = 1.0
    let GalerkinSolution = GalerkinMethod equation a b n
    let (networkSolution, errors, nIterations) = solveWithNetwork equation a b n 1e-6 7
    let n1 = networkSolution.Length
    let x = buildNetwork a b (n1 - 1)
    let expected = List.init n1 (fun i -> u x[i])
    let actualGalerkin = List.init n1 (fun i -> GalerkinSolution x[i])
    printfn "%A" actualGalerkin
    drawPlot3 x expected (networkSolution |> Vector.toArray |> Array.toList) actualGalerkin

test1 ()
