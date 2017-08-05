/**
 * 
 */
package util;


/**
 * Represents an operation that accepts no arguments and returns no result.
 * Unlike most other functional interfaces, BiConsumer is expected to operate
 * via side-effects.
 * 
 * This is a functional interface whose functional method is execute().
 * 
 * I didn't want to use Runnable because of its association with Threads, so I
 * just made my own functional interface.
 * 
 * @author jkunimune
 */
public interface Procedure {
	public static final Procedure NONE = ()->{};

	public abstract void execute();
	
	public static Procedure concat(Procedure...procedures) {
		return () -> {
				for (Procedure p: procedures)
					p.execute();
			};
	}
}
