
class Runnable
  def cond
    if @cond == nil
      @cond = self.new_cond
    end
    return @cond
  end

  def start()
    @thread_obj = Thread.new(self){|r| r.run()}
  end
end
