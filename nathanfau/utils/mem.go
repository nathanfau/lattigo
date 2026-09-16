package utils

// Memory helpers for instrumenting a pipeline: what the Go heap keeps alive, what the process holds,
// and the most it ever held. Nothing CKKS-specific.

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"runtime"
	"runtime/debug"
	"strconv"
	"strings"
)

// MemSnapshot is where the memory stands at one point of a run.
type MemSnapshot struct {
	Heap    uint64 // Go heap in use: live objects plus the garbage not collected yet
	HeapSys uint64 // heap the runtime holds from the OS, in use or not
	RSS     uint64 // resident memory of the process; 0 where /proc is missing
	PeakRSS uint64 // highest RSS so far (VmHWM): what the machine actually had to provide
}

// Mem reads the current snapshot without collecting: cheap enough for a timed section, but its
// Heap still carries the garbage of the last operations. LiveHeap is the collected figure.
func Mem() MemSnapshot {
	var ms runtime.MemStats
	runtime.ReadMemStats(&ms)
	s := MemSnapshot{Heap: ms.HeapAlloc, HeapSys: ms.HeapSys}
	s.RSS, s.PeakRSS = procStatus()
	return s
}

// LiveHeap collects first, so it returns what the program keeps alive. A cycle is cheap on
// lattigo's heap, whose bulk is pointer-free []uint64 the collector does not scan, but it is still a
// cycle: keep it out of timed sections.
func LiveHeap() uint64 {
	runtime.GC()
	var ms runtime.MemStats
	runtime.ReadMemStats(&ms)
	return ms.HeapAlloc
}

func (s MemSnapshot) String() string {
	if s.RSS == 0 {
		return fmt.Sprintf("heap %s (held %s)", Bytes(s.Heap), Bytes(s.HeapSys))
	}
	return fmt.Sprintf("heap %s (held %s), RSS %s, peak %s", Bytes(s.Heap), Bytes(s.HeapSys), Bytes(s.RSS), Bytes(s.PeakRSS))
}

// GCSettings names the knobs that decide how far the heap may grow past what is alive: with
// GOGC=100 and no GOMEMLIMIT, the runtime lets the heap reach twice the live heap before collecting.
func GCSettings() string {
	gogc := os.Getenv("GOGC")
	if gogc == "" {
		gogc = "100"
	}
	limit := "none"
	if l := debug.SetMemoryLimit(-1); l != math.MaxInt64 {
		limit = Bytes(uint64(l))
	}
	return fmt.Sprintf("GOGC=%s, GOMEMLIMIT=%s, GOMAXPROCS=%d", gogc, limit, runtime.GOMAXPROCS(0))
}

// Bytes formats a size in Ko, Mo, Go, To, one unit being 1024 of the one below: 1 Go = 1024 Mo =
// 2^30 octets. Every size the pipeline prints goes through here, so they all read the same way.
func Bytes(n uint64) string {
	const unit = 1 << 10
	if n < unit {
		return fmt.Sprintf("%d o", n)
	}
	div, exp := uint64(unit), 0
	for m := n / unit; m >= unit; m /= unit {
		div *= unit
		exp++
	}
	return fmt.Sprintf("%.2f %co", float64(n)/float64(div), "KMGTPE"[exp])
}

// BytesDelta formats after - before, signed.
func BytesDelta(after, before uint64) string {
	if after >= before {
		return "+" + Bytes(after-before)
	}
	return "-" + Bytes(before-after)
}

// ResetPeakRSS brings the peak RSS back to the current RSS, so a run of several configurations can
// read each one's own peak. Linux only (/proc/self/clear_refs, value 5); a no-op elsewhere.
func ResetPeakRSS() {
	_ = os.WriteFile("/proc/self/clear_refs", []byte("5"), 0)
}

// procStatus reads VmRSS and VmHWM from /proc/self/status, in bytes; zeros where there is no /proc.
func procStatus() (rss, peak uint64) {
	f, err := os.Open("/proc/self/status")
	if err != nil {
		return 0, 0
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		key, val, ok := strings.Cut(sc.Text(), ":")
		if !ok {
			continue
		}
		var dst *uint64
		switch key {
		case "VmRSS":
			dst = &rss
		case "VmHWM":
			dst = &peak
		default:
			continue
		}
		if fields := strings.Fields(val); len(fields) > 0 { // "123456 kB"
			if kb, err := strconv.ParseUint(fields[0], 10, 64); err == nil {
				*dst = kb << 10
			}
		}
	}
	return rss, peak
}
