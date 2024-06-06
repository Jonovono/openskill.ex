defmodule OpenskillStatistics do
  @moduledoc false

  @sqrt_2_pi :math.sqrt(2 * :math.pi())
  @inv_sqrt_2 :math.sqrt(0.5)

  def phi_major(x) do
    0.5 * (1.0 + :math.erf(x / :math.sqrt(2.0)))
  end

  def phi_major_inverse(x) do
    :math.sqrt(2) * :math.erfinv(2 * x - 1)
  end

  def phi_minor(x) do
    :math.exp(-0.5 * x * x) / @sqrt_2_pi
  end

  def v(x, t) do
    xt = x - t
    denom = phi_major(xt)
    if denom < :math.pow(2, -52) do
      -xt
    else
      phi_minor(xt) / denom
    end
  end

  @spec w(number(), number()) :: number()
  def w(x, t) do
    xt = x - t
    denom = phi_major(xt)
    if denom < :math.pow(2, -52) do
      if x < 0 do
        1
      else
        0
      end
    else
      v(x, t) * (v(x, t) + xt)
    end
  end

  def vt(x, t) do
    xx = abs(x)
    b = phi_major(t - xx) - phi_major(-t - xx)
    if b < 1.0e-5 do
      if x < 0 do
        -x - t
      else
        -x + t
      end
    else
      a = phi_minor(-t - xx) - phi_minor(t - xx)
      if x < 0 do
        -a / b
      else
        a / b
      end
    end
  end

  def wt(x, t) do
    xx = abs(x)
    b = phi_major(t - xx) - phi_major(-t - xx)
    if b < :math.pow(2, -52) do
      1.0
    else
      ((t - xx) * phi_minor(t - xx) + (t + xx) * phi_minor(-t - xx)) / b + vt(x, t) * vt(x, t)
    end
  end
end
