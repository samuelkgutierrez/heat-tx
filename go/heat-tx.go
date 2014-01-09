// A simple 2D heat transfer simulation in Go by Samuel K. Gutierrez

// To Profile:
// Build
// ./heat-tx -cpuprofile=heat-tx.prof
// go tool pprof ./heat-tx ./heat-tx.prof

// To Plot (gnuplot):
// plot './heat-img.dat' matrix with image

package main

import (
    "flag"
    "strconv"
    "fmt"
    "bufio"
    "os"
    "log"
    "runtime/pprof"
)

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")

// Application constants
const (
    // Application Name
    AppName string = "go-heat-tx"
    // Application version string
    AppVerStr string = "0.1"
    // Mesh size (x and y)
    N uint64 = 512
    // Max simulation time
    TMax uint64 = 1024
    // Some constant
    K float64 = 0.4
    // Thermal conductivity
    ThermCond float64 = 0.6
)

// 2D mesh
type Mesh struct {
    nx, ny uint64 // mesh size in x and y
    cells [][]float64 // mesh cells
}

type SimParams struct {
    // Thermal conductivity
    c float64
    deltaS float64
    // time interval
    deltaT float64
    // Max sim iterations
    tMax uint64
}

type HeatTxSim struct {
    // The meshes
    newMesh, oldMesh *Mesh
    // Simulation parameters
    params *SimParams
}

// NewMesh returns an empty mesh of the specified width and height
func NewMesh(x, y uint64) *Mesh {
    cells := make([][]float64, x)
    for i := range cells {
        cells[i] = make([]float64, y)
    }
    // **remember** unlike C, we can return the address of a local variable. 
    // in fact, this returns a fresh instance each time the following code is
    // evaluated - w00t.
    return &Mesh{nx: x, ny: y, cells: cells}
}

// NewSimParams returns a new set of initialized simulation parameters based on
// the provided input.
func NewSimParams(nx uint64, thermCond float64, tMax uint64) *SimParams {
    fmt.Println("o initializing simulation parameters...")
    ds := 1.0 / (float64(nx) + 1.0)
    dt := (ds * ds) / (4.0 * thermCond)
    sp := &SimParams{c: thermCond, deltaS: ds, deltaT: dt, tMax: tMax}
    fmt.Print(sp)
    return sp
}

// Nice SimParams printing
func (p *SimParams) String() string {
    pStr := ""
    pStr += fmt.Sprintf(". max_t: %d\n", p.tMax)
    pStr += fmt.Sprintf(". c: %f\n", p.c)
    pStr += fmt.Sprintf(". delta_s: %f\n", p.deltaS)
    pStr += fmt.Sprintf(". delta_t: %f\n", p.deltaT)
    pStr += "\n"
    return pStr
}

// Nice Mesh String representation
func (m *Mesh) String() string {
    mStr := ""
    for i := range m.cells {
        for j := range m.cells[i] {
            mStr += strconv.FormatFloat(m.cells[i][j], 'e', 1, 64) + " "
        }
        mStr += "\n"
    }
    return mStr
}

func NewHeatTxSim(x, y uint64,  thermCond float64, tMax uint64) *HeatTxSim {
    return &HeatTxSim{params: NewSimParams(x, thermCond, tMax),
                      newMesh: NewMesh(x, y), oldMesh: NewMesh(x, y)}
}

func (m *Mesh) SetInitConds() {
    x0 := m.nx / 2
    y0 := m.ny / 2
    x  := m.nx / 4
    y  := uint64(0)
    radiusErr := int64(1 - x)

    for x >= y {
        m.cells[ x + x0][ y + y0] = K
        m.cells[ x + x0][ y + y0] = K * .50
        m.cells[ y + x0][ x + y0] = K * .60
        m.cells[-x + x0][ y + y0] = K * .70
        m.cells[-y + x0][ x + y0] = K * .80
        m.cells[-x + x0][-y + y0] = K * .70
        m.cells[-y + x0][-x + y0] = K * .60
        m.cells[ x + x0][-y + y0] = K * .50
        m.cells[ y + x0][-x + y0] = K
        y++
        if radiusErr < 0 {
            radiusErr += int64(2 * y + 1)
        } else {
            x--
            radiusErr += int64(2 * (y - x + 1))
        }
    }
}

// Runs the simulation
func (s *HeatTxSim) Run() {
    fmt.Println("o starting simulation...")
    nx := len(s.oldMesh.cells) - 1
    ny := len(s.oldMesh.cells[0]) - 1
    ds2 := s.params.deltaS*s.params.deltaS
    cdtods2 := s.params.c * s.params.deltaT / ds2
    // Stash the mesh pointers
    newMesh := s.newMesh
    oldMesh := s.oldMesh

    for t := uint64(0); t < s.params.tMax; t++ {
        if t % 100 == 0 {
            fmt.Println(". starting iteration", t, "of", s.params.tMax)
        }
        for i := 1; i < nx; i++ {
            for j := 1; j < ny; j++ {
                newMesh.cells[i][j] =
                    oldMesh.cells[i][j] +
                    (cdtods2 *
                     (oldMesh.cells[i + 1][j] +
                      oldMesh.cells[i - 1][j] -
                      4.0 * oldMesh.cells[i][j] +
                      oldMesh.cells[i][j + 1] +
                      oldMesh.cells[i][j - 1]))
            }
        }
        // swap old and new - this is just a pointer swap
        oldMesh, newMesh = newMesh, oldMesh
        // Constant heat source
        oldMesh.SetInitConds()
    }
}

// Dumps to a text file with the current newMesh cell values
func (s *HeatTxSim) Dump() error {
    dumpFile, err := os.Create("heat-img.dat")
    if err != nil { return err }
    defer func() {
        if err := dumpFile.Close(); err != nil {
            panic(err)
        }
    }()
    // Create a new buffered writer
    w := bufio.NewWriter(dumpFile)
    for i := range s.newMesh.cells {
        for j := range s.newMesh.cells[0] {
            fmt.Fprintf(w, "%f", s.newMesh.cells[i][j])
            if j != len(s.newMesh.cells) - 1 {
                fmt.Fprintf(w, " ")
            }
        }
        fmt.Fprintln(w)
    }
    return w.Flush()
}

func main() {
    flag.Parse()
    if *cpuprofile != "" {
        f, err := os.Create(*cpuprofile)
        if err != nil {
            log.Fatal(err)
        }
        pprof.StartCPUProfile(f)
        defer pprof.StopCPUProfile()
    }
    fmt.Println("o", AppName, AppVerStr)
    sim := NewHeatTxSim(N, N, ThermCond, TMax)
    sim.oldMesh.SetInitConds()
    sim.Run()
    err := sim.Dump()
    if (err != nil) { panic(err) }
}
