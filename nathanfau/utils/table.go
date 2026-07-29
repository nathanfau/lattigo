package utils

import "strings"

// TableGroup is a set of sub-columns sharing one merged header in a BoxTable.
type TableGroup struct {
	Name string
	Cols []string // sub-column headers, left to right
}

// BoxTable renders an aligned table with Unicode box-drawing borders, in the style of the CKKS
// precision-stats table. The first column is an index column headed by lead; the other columns are
// organised into groups, each carrying a header that spans its sub-columns. leads[i] is row i's
// index value and rows[i] holds that row's sub-column values flattened left-to-right across every
// group, so len(rows[i]) should equal the total number of sub-columns. Column widths are computed
// to fit their content; headers are centered and values right-aligned. Cell contents are assumed
// to be single-width ASCII. The returned string ends with a newline.
//
//	BoxTable("start",
//	    []TableGroup{{"Cleaning", []string{"prec", "time"}}, {"Smoother", []string{"prec", "time"}}},
//	    []string{"2.00", "4.00"},
//	    [][]string{{"2.43", "20ms", "2.72", "30ms"}, {"6.42", "20ms", "8.69", "31ms"}})
func BoxTable(lead string, groups []TableGroup, leads []string, rows [][]string) string {
	// Flatten the sub-headers and record each group's span into the flattened columns.
	var subs []string
	starts := make([]int, len(groups))
	for g, grp := range groups {
		starts[g] = len(subs)
		subs = append(subs, grp.Cols...)
	}

	// Content widths: the longest of the header and every value in each column.
	lw := len(lead)
	for _, v := range leads {
		lw = max(lw, len(v))
	}
	cw := make([]int, len(subs))
	for j, h := range subs {
		cw[j] = len(h)
	}
	for _, row := range rows {
		for j := 0; j < len(subs) && j < len(row); j++ {
			cw[j] = max(cw[j], len(row[j]))
		}
	}

	// Field widths add one space of padding on each side.
	flw := lw + 2
	fw := make([]int, len(subs))
	for j := range fw {
		fw[j] = cw[j] + 2
	}

	// groupWidth is the display width a group header spans: its sub-columns plus their inner
	// dividers. Widen the last sub-column if the group name does not fit.
	groupWidth := func(g int) int {
		w := len(groups[g].Cols) - 1
		for k := range groups[g].Cols {
			w += fw[starts[g]+k]
		}
		return w
	}
	for g := range groups {
		if need := len(groups[g].Name) + 2; need > groupWidth(g) {
			fw[starts[g]+len(groups[g].Cols)-1] += need - groupWidth(g)
		}
	}

	fill := func(n int) string { return strings.Repeat("─", n) }
	var b strings.Builder

	// Top border.
	b.WriteString("┌")
	b.WriteString(fill(flw))
	for g := range groups {
		b.WriteString("┬")
		b.WriteString(fill(groupWidth(g)))
	}
	b.WriteString("┐\n")

	// Merged header row: index header + group names.
	b.WriteString("│")
	b.WriteString(center(lead, flw))
	for g := range groups {
		b.WriteString("│")
		b.WriteString(center(groups[g].Name, groupWidth(g)))
	}
	b.WriteString("│\n")

	// Divider: the index column stays merged (blank), the groups split into their sub-columns.
	b.WriteString("│")
	b.WriteString(strings.Repeat(" ", flw))
	b.WriteString("├")
	for g := range groups {
		if g > 0 {
			b.WriteString("┼")
		}
		for k := range groups[g].Cols {
			if k > 0 {
				b.WriteString("┬")
			}
			b.WriteString(fill(fw[starts[g]+k]))
		}
	}
	b.WriteString("┤\n")

	// Sub-header row.
	b.WriteString("│")
	b.WriteString(strings.Repeat(" ", flw))
	for j, h := range subs {
		b.WriteString("│")
		b.WriteString(center(h, fw[j]))
	}
	b.WriteString("│\n")

	// Mid border: every column is now separate.
	b.WriteString("├")
	b.WriteString(fill(flw))
	for j := range subs {
		b.WriteString("┼")
		b.WriteString(fill(fw[j]))
	}
	b.WriteString("┤\n")

	// Data rows.
	for i, row := range rows {
		lv := ""
		if i < len(leads) {
			lv = leads[i]
		}
		b.WriteString("│")
		b.WriteString(rjustCell(lv, flw))
		for j := range subs {
			v := ""
			if j < len(row) {
				v = row[j]
			}
			b.WriteString("│")
			b.WriteString(rjustCell(v, fw[j]))
		}
		b.WriteString("│\n")
	}

	// Bottom border.
	b.WriteString("└")
	b.WriteString(fill(flw))
	for j := range subs {
		b.WriteString("┴")
		b.WriteString(fill(fw[j]))
	}
	b.WriteString("┘\n")

	return b.String()
}

// center pads s with spaces to width w (extra space goes to the right).
func center(s string, w int) string {
	if len(s) >= w {
		return s[:w]
	}
	left := (w - len(s)) / 2
	return strings.Repeat(" ", left) + s + strings.Repeat(" ", w-len(s)-left)
}

// rjustCell right-aligns s in a field of width w, keeping at least one space of padding each side.
func rjustCell(s string, w int) string {
	if len(s) > w-2 {
		if len(s) > w {
			return s[:w]
		}
		return s
	}
	return strings.Repeat(" ", w-1-len(s)) + s + " "
}
