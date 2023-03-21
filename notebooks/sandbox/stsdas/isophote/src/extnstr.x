# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
# EXTNSTR -- Extracts n-th substring from dictionary list. The same
# format is used as in strdic.

procedure extnstr (str, n, substr)

char	str[], substr[]
int	n

int	i, j, ip1, ip2, strnidx()

begin
	ip1 = strnidx (str[1], n,   str)
	ip2 = strnidx (str[1], n+1, str)
	j = 1
	do i = ip1+1, ip2-1 {
	    substr[j] = str[i]
	    j = j + 1
	}
	substr[j] = EOS
end


# STRNIDX -- Return the index of the n-th occurrence of a character in a
# string.

int procedure strnidx (ch, n, str)

char	ch, str[ARB]
int	n, i, ip

begin
	i = 0
	do ip = 1, ARB {
	    if (str[ip] == EOS)
	        return (ip)
	    else if (str[ip] == ch) {
	        i = i + 1
	        if (i >= n)
	            return (ip)
	    }
	}
end
