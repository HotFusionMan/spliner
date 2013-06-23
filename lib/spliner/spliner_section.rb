require 'matrix'
require 'spliner/spliner_section'

module Spliner

  # Spliner::SplinerSection is only used via Spliner::Spliner
  #
  # As the spline algorithm does not handle duplicate X values well, the 
  # curve is split into two non continuous parts where duplicate X values 
  # appear. Each such part is represented by a SplinerSection
  class SplinerSection
    unless RUBY_VERSION >= "1.9.3"
      # Needed for the LUP decomposition code that was copy-pasted in from Ruby 2.0.0-p195's Matrix.
      include Matrix::ConversionHelper
    end

    attr_reader :k, :x, :y

    def initialize(x, y)
      @x, @y = x, y
      @x_pairs = @x.each_cons(2).map {|pair| pair.first..pair.last }
      check_points_increasing
      calculate_a_k
    end

    def range
      @x.first..@x.last
    end

    def calculate_a_k
      if @x.size > 1
        inv_diff = @x.each_cons(2).map {|x1, x2| 1.0 / (x2 - x1) }
        a_diag = 2.0 * Matrix::diagonal(*vector_helper(inv_diff))
        a_non_diag = Matrix::build(@x.size) do |row, col|
          if row == col + 1
            inv_diff[col]
          elsif col == row + 1
            inv_diff[row]
          else
            0.0
          end
        end

        a = a_diag + a_non_diag

        tmp = @x.zip(@y).each_cons(2).map do |p1, p2|
          x1, y1 = p1
          x2, y2 = p2
          delta_x = (x2 - x1)
          3.0 * (y2 - y1) / (delta_x * delta_x)
        end
        b = vector_helper(tmp)

        if RUBY_VERSION >= "1.9.3"
          @k = a.lup_decomposition.solve b
        else
          # Although a.inverse * b actually passes the specs to within 2 * epsilon precision,
          # LUP decomposition is a better algorithm in general, so copy-paste it here from Ruby 2.0.0's Matrix library.

          @lu = a.to_a
          @row_size = a.row_size
          @col_size = a.column_size
          @pivots = Array.new(@row_size)
          @row_size.times do |i|
            @pivots[i] = i
          end
          @pivot_sign = 1
          lu_col_j = Array.new(@row_size)

          # Outer loop.

          @col_size.times do |j|

            # Make a copy of the j-th column to localize references.

            @row_size.times do |i|
              lu_col_j[i] = @lu[i][j]
            end

            # Apply previous transformations.

            @row_size.times do |i|
              lu_row_i = @lu[i]

              # Most of the time is spent in the following dot product.

              kmax = [i, j].min
              s = 0
              kmax.times do |k|
                s += lu_row_i[k]*lu_col_j[k]
              end

              lu_row_i[j] = lu_col_j[i] -= s
            end

            # Find pivot and exchange if necessary.

            p = j
            (j+1).upto(@row_size-1) do |i|
              if (lu_col_j[i].abs > lu_col_j[p].abs)
                p = i
              end
            end
            if (p != j)
              @col_size.times do |k|
                t = @lu[p][k]; @lu[p][k] = @lu[j][k]; @lu[j][k] = t
              end
              k = @pivots[p]; @pivots[p] = @pivots[j]; @pivots[j] = k
              @pivot_sign = -@pivot_sign
            end

            # Compute multipliers.

            if (j < @row_size && @lu[j][j] != 0)
              (j+1).upto(@row_size-1) do |i|
                @lu[i][j] = @lu[i][j].quo(@lu[j][j])
              end
            end
          end


          # And now for the backsolve to get @k:
          b = convert_to_array(b)
          if (b.size != @row_size)
            Matrix.Raise Matrix::ErrDimensionMismatch
          end

          # Copy right hand side with pivoting
          m = b.values_at(*@pivots)

          # Solve L*Y = P*b
          @col_size.times do |k|
            (k+1).upto(@col_size-1) do |i|
              m[i] -= m[k]*@lu[i][k]
            end
          end
          # Solve U*m = Y
          (@col_size-1).downto(0) do |k|
            m[k] = m[k].quo(@lu[k][k])
            k.times do |i|
              m[i] -= m[k]*@lu[i][k]
            end
          end
          @k = Vector.elements(m, false)
        end
      else
        @k = Vector[0.0]
      end
    end
    private :calculate_a_k

    # returns an interpolated value
    def get(v)
      i = @x_pairs.find_index {|pair| pair.cover? v }
      if i
        x_pair = @x_pairs[i]
        x_min = x_pair.min
        dx = x_pair.max - x_min
        y_max = @y[i + 1]
        y_min = @y[i]
        dy = y_max - y_min
        t = (v - x_min) / dx
        a =      @k[i] * dx - dy
        b = -@k[i + 1] * dx + dy
        one_minus_t = 1 - t
        t * y_max + one_minus_t * ( y_min + t * ( a * one_minus_t + b * t ) )
      elsif @x.size == 1 && @x.first == v
        @y.first
      else
        nil
      end
    end

    # for a vector [a, b, c] returns [a, a + b, b + c, c]
    def vector_helper(a)
      Vector[*([0.0] + a)] + Vector[*(a + [0.0])]
    end
    private :vector_helper

    def check_points_increasing
      @x.each_cons(2) do |x1, x2|
        raise 'Key point\'s X values should be in increasing order' unless x2 > x1
      end
    end
    private :check_points_increasing
  end
end

