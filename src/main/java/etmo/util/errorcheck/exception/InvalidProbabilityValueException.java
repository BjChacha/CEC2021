package etmo.util.errorcheck.exception;

@SuppressWarnings("serial")
public class InvalidProbabilityValueException extends RuntimeException {
  public InvalidProbabilityValueException(double value) {
    super("The parameter " + value + " is not a valid probability value") ;
  }
}
