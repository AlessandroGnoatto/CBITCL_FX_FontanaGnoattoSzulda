/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 11.10.2013
 */

package org.calibrationframework.modelling;

import java.util.Map;

/**
 * Interface implemented by all financial product which may be valued by a model.
 *
 * @author Christian Fries
 */
public interface ProductInterface {

	/**
	 * Return the valuation of the product using the given model.
	 *
	 * Implement this method using a checked cast of the model to a derived model for which the product
	 * provides a valuation algorithm. Example: an interest rate product requires that the passed model
	 * object implements the interface of an interest rate model. Since there is no polymorphism on
	 * arguments (see Double Dynamic Dispatch), we reply on a checked cast.
	 *
	 * @param model The model under which the product is valued.
	 * @return Object containing the value of the product using the given model.
	 */
	Object getValue(ModelInterface model);

	/**
	 * Return the valuation of the product using the given model.
	 *
	 * Implement this method using a checked cast of the model to a derived model for which the product
	 * provides a valuation algorithm. Example: an interest rate product requires that the passed model
	 * object implements the interface of an interest rate model. Since there is no polymorphism on
	 * arguments (see Double Dynamic Dispatch), we reply on a checked cast.
	 *
	 * @param model The model under which the product is valued.
	 * @return Map containing the value of the product using the given model.
	 */
	default Map<String, Object> getValues(ModelInterface model) throws UnsupportedOperationException { throw new UnsupportedOperationException(); }
}
