base10exp(v::AbstractFloat) = v==0 ? 0 : round(Int, floor(log10(abs(v))))

function parts10(v::AbstractFloat, exp_range)
  e = base10exp(v)
  if e in exp_range
    e = 0
  end

  m = v / exp10(e)
  return (m, e)
end

function format(v)
  string(v)
end

# I don't really trust Base.trunc to do base-10 operations? but I should.
function trunc_str(v::AbstractString, at)
  if M[1] == '-' at += 1 end # handle minus sign
  # round last digit upward (maybe should be away from zero)
  last_digit = parse(Int, M[at]) + (parse(Int, M[at+1]) >= 5)
  return M[1:(at-1)] * string(last_digit)
end

# try to print float64s the way that fortran does
function format(v::AbstractFloat)
  (m, e) = parts10(v, -1:18)
  #M = @sprintf("%.16f", m)
  #M = M[1]=='-' ?  M[1:19] : M[1:18] # handle minus sign
  if abs(m) >= 1 || m == 0
    M = @sprintf("%.16f", m)
  else # first digit is zero, not significant.
    M = @sprintf("%.17f", m)
  end
  es = e >= 0 ? "+" : "-"
  E = "E" * es * @sprintf("%03d", abs(e))
  E0 = e != 0 ? E : ""
  return M*E0
end

function format(av::Array)
  join(map(format, av), " ")
end

function format(args...)
  join(map(format, args), " ")
end

function fortranprint(args...)
  print(" ")
  println(format(args...))
end
