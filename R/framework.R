

set_level_mean = function(s, l) {
	mean(s[l], na.rm = TRUE)
}

set_level_ks = function(s, l) {

	od = order(s, decreasing = TRUE)
	s = s[od]
	l = l[od]

	n = length(s)

    s_set = numeric(n)
    s_set[l] = abs(s[l])

    l_other = !l

    f1 = cumsum(s_set)/sum(s_set)
    f2 = cumsum(l_other)/sum(l_other)
    m1 = max(f1 - f2)
    m2 = min(f1 - f2)

    if(m1 >= 0 && m2 >= 0) {
        es = max(m1, m2)
    } else if(m1 <= 0 && m2 <= 0) {
        es = min(m1, m2)
    } else {
        es = m1 + m2
    }  

    es
}

